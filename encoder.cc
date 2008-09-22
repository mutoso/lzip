/*  Lzip - A LZMA file compressor
    Copyright (C) 2008 Antonio Diaz Diaz.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#define _FILE_OFFSET_BITS 64

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <stdint.h>

#include "lzip.h"
#include "encoder.h"


const Prob_prices prob_prices;


bool Matchfinder::read_block() throw()
  {
  const int size = buffer_size - stream_pos;
  if( size == 0 ) return true;
  const int rd = readblock( _ides, (char *)buffer + stream_pos, size );
  if( rd != size && errno ) return false;

  stream_pos += rd;
  if( rd == size ) pos_limit = stream_pos - after_size;	// last safe position
  else { pos_limit = stream_pos; at_stream_end = true; }
  return true;
  }


bool Matchfinder::move_pos() throw()
  {
  update_crc( _crc, buffer[pos] );
  if( ++cyclic_pos >= dictionary_size ) cyclic_pos = 0;
  if( ++pos >= pos_limit )
    {
    if( !at_stream_end )
      {
      const int offset = pos - dictionary_size - max_num_trials;
      const int size = stream_pos - offset;
      std::memmove( buffer, buffer + offset, size );
      partial_file_pos += offset;
      pos -= offset;
      stream_pos -= offset;
      pos_limit -= offset;
      for( int i = 0; i < num_prev_positions; ++i )
        if( prev_positions[i] >= 0 ) prev_positions[i] -= offset;
      for( int i = 0; i < dictionary_size; ++i )
        if( prev_pos_tree[i] >= 0 ) prev_pos_tree[i] -= offset;
      return read_block();
      }
    else if( pos > pos_limit ) { pos = pos_limit; return false; }
    }
  return true;
  }


int Matchfinder::longest_match_len( int * const distances ) throw()
  {
  int len_limit = match_len_limit;
  if( len_limit > available_bytes() )
    {
    len_limit = available_bytes();
    if( len_limit < min_match_len ) return 0;
    }

  const uint8_t * const data = ptr_to_current_pos();
  const int key = ( (int)data[0] << 8 ) + data[1];

  int newpos = prev_positions[key];
  prev_positions[key] = pos;

  const int min_pos = (pos > dictionary_size) ? (pos - dictionary_size) : 0;
  int maxlen = min_match_len - 1;
  int * ptr0 = prev_pos_tree + ( ( cyclic_pos << 1 ) & dictionary_mask );
  int * ptr1 = ptr0 + 1;

  for( int count = 256;; )
    {
    if( newpos < min_pos ) { *ptr0 = *ptr1 = -1; break; }
    const uint8_t * const newdata = buffer + newpos;
    if( newdata[0] != data[0] || newdata[1] != data[1] ) break;
    int len = 2;
    while( len < len_limit && newdata[len] == data[len] ) ++len;

    const int delta = pos - newpos;
    if( distances ) while( maxlen < len ) distances[++maxlen] = delta - 1;

    int * const newptr = prev_pos_tree + ( ( ( cyclic_pos - delta ) << 1 ) & dictionary_mask );

    if( len < len_limit )
      {
      if( newdata[len] < data[len] )
        { *ptr0 = newpos; ptr0 = newptr + 1; newpos = *ptr0; }
      else
        { *ptr1 = newpos; ptr1 = newptr;     newpos = *ptr1; }
      }
    else
      { *ptr0 = newptr[0]; *ptr1 = newptr[1]; break; }
    if( --count <= 0 ) break;
    }
  return maxlen;
  }


int Len_encoder::calculate_price( const int symbol, const int pos_state ) const throw()
  {
  if( symbol < len_low_symbols )
    return price0( choice1 ) + price_symbol( bm_low[pos_state], symbol, len_low_bits );
  int price = price1( choice1 );
  if( symbol < len_low_symbols + len_mid_symbols )
    {
    price += price0( choice2 );
    price += price_symbol( bm_mid[pos_state], symbol - len_low_symbols, len_mid_bits );
    }
  else
    {
    price += price1( choice2 );
    price += price_symbol( bm_high, symbol - len_low_symbols - len_mid_symbols, len_high_bits );
    }
  return price;
  }


void Len_encoder::encode( Range_encoder & range_encoder, const int symbol,
                          const int pos_state )
  {
  if( symbol < len_low_symbols )
    {
    range_encoder.encode_bit( choice1, 0 );
    range_encoder.encode_tree( bm_low[pos_state], symbol, len_low_bits );
    }
  else
    {
    range_encoder.encode_bit( choice1, 1 );
    if( symbol < len_low_symbols + len_mid_symbols )
      {
      range_encoder.encode_bit( choice2, 0 );
      range_encoder.encode_tree( bm_mid[pos_state], symbol - len_low_symbols, len_mid_bits );
      }
    else
      {
      range_encoder.encode_bit( choice2, 1 );
      range_encoder.encode_tree( bm_high, symbol - len_low_symbols - len_mid_symbols, len_high_bits );
      }
    }
  if( --counters[pos_state] <= 0 ) update_prices( pos_state );
  }


void LZ_encoder::fill_align_prices() throw()
  {
  for( int i = 0; i < dis_align_size; ++i )
    align_prices[i] = price_symbol_reversed( bm_align, i, dis_align_bits );
  align_price_count = dis_align_size;
  }


void LZ_encoder::fill_distance_prices() throw()
  {
  for( int dis_state = 0; dis_state < max_dis_states; ++dis_state )
    {
    int slot = 0;
    for( ; slot < end_dis_model && slot < num_dis_slots; ++slot )
      dis_slot_prices[dis_state][slot] =
        price_symbol( bm_dis_slot[dis_state], slot, dis_slot_bits );
    for( ; slot < num_dis_slots; ++slot )
      dis_slot_prices[dis_state][slot] =
        price_symbol( bm_dis_slot[dis_state], slot, dis_slot_bits ) +
        (((( slot >> 1 ) - 1 ) - dis_align_bits ) << price_shift );

    int dis = 0;
    for( ; dis < start_dis_model; ++dis )
      dis_prices[dis_state][dis] = dis_slot_prices[dis_state][dis];
    for( ; dis < modeled_distances; ++dis )
      {
      const int dis_slot = get_dis_slot( dis );
      const int direct_bits = ( dis_slot >> 1 ) - 1;
      const int base = ( 2 | ( dis_slot & 1 ) ) << direct_bits;
      dis_prices[dis_state][dis] = dis_slot_prices[dis_state][dis_slot] +
          price_symbol_reversed( bm_dis + base - dis_slot, dis - base, direct_bits );
      }
    }
  }


// Return value: ( dis == -1 ) && ( len == 1 ) means literal
int LZ_encoder::best_pair_sequence( const int reps[num_rep_distances],
                                    const State & state )
  {
  int main_len;
  if( longest_match_found > 0 )		// from previous call
    {
    main_len = longest_match_found;
    longest_match_found = 0;
    }
  else main_len = read_match_distances();

  int replens[num_rep_distances];
  int rep_index = 0;
  for( int i = 0; i < num_rep_distances; ++i )
    {
    replens[i] = matchfinder.true_match_len( 0, reps[i] + 1, max_match_len );
    if( replens[i] > replens[rep_index] ) rep_index = i;
    }
  if( replens[rep_index] >= match_len_limit )
    {
    trials[0].dis = rep_index;
    trials[0].price = replens[rep_index];
    if( !move_pos( replens[rep_index] ) ) return 0;
    return replens[rep_index];
    }

  if( main_len >= match_len_limit )
    {
    trials[0].dis = match_distances[match_len_limit] + num_rep_distances;
    trials[0].price = main_len;
    if( !move_pos( main_len ) ) return 0;
    return main_len;
    }

  trials[0].state = state;
  for( int i = 0; i < num_rep_distances; ++i ) trials[0].reps[i] = reps[i];

  const uint8_t cur_byte = matchfinder[0];
  const uint8_t match_byte = matchfinder[-reps[0]-1];
  unsigned int position = matchfinder.file_position();
  const int pos_state = position & pos_state_mask;

  trials[1].dis = -1;
  trials[1].prev_index = 0;
  trials[1].price = price0( bm_match[state()][pos_state] );
  if( state.is_char() )
    trials[1].price += literal_encoder.price_symbol( matchfinder[-1], cur_byte );
  else
    trials[1].price += literal_encoder.price_matched( matchfinder[-1], cur_byte, match_byte );

  const int match_price = price1( bm_match[state()][pos_state] );
  const int rep_match_price = match_price + price1( bm_rep[state()] );

  if( match_byte == cur_byte )
    trials[1].update( 0, 0, rep_match_price + price_rep_len1( state, pos_state ) );

  if( main_len < min_match_len )
    {
    trials[0].dis = trials[1].dis;
    trials[0].price = 1;
    if( !matchfinder.move_pos() ) return 0;
    return 1;
    }

  {
  const int normal_match_price = match_price + price0( bm_rep[state()] );
  int len = min_match_len;
  if( main_len <= replens[rep_index] )
    {
    main_len = replens[rep_index];
    for( ; len <= main_len; ++len ) trials[len].price = infinite_price;
    }
  else for( ; len <= main_len; ++len )
    {
    trials[len].dis = match_distances[len] + num_rep_distances;
    trials[len].prev_index = 0;
    trials[len].price = normal_match_price +
                        price_pair( match_distances[len], len, pos_state );
    }
  }

  for( int rep = 0; rep < num_rep_distances; ++rep )
    for( int len = min_match_len; len <= replens[rep]; ++len )
      trials[len].update( rep, 0, rep_match_price +
                                  price_rep( rep, len, state, pos_state ) );

  int cur = 0;
  int num_trials = main_len;
  if( !matchfinder.move_pos() ) return 0;

  while( true )
    {
    if( ++cur >= num_trials )
      {
      backward( cur );
      return cur;
      }
    const int newlen = read_match_distances();
    if( newlen >= match_len_limit )
      {
      longest_match_found = newlen;
      backward( cur );
      return cur;
      }

    Trial & cur_trial = trials[cur];
    const int prev_index = cur_trial.prev_index;

    cur_trial.state = trials[prev_index].state;

    for( int i = 0; i < num_rep_distances; ++i )
      cur_trial.reps[i] = trials[prev_index].reps[i];
    if( prev_index == cur - 1 )
      {
      if( cur_trial.dis == 0 ) cur_trial.state.set_short_rep();
      else cur_trial.state.set_char();
      }
    else
      {
      if( cur_trial.dis < num_rep_distances ) cur_trial.state.set_rep();
      else cur_trial.state.set_match();
      mtf_reps( cur_trial.dis, cur_trial.reps );
      }

    const uint8_t cur_byte = matchfinder[0];
    const uint8_t match_byte = matchfinder[-cur_trial.reps[0]-1];
    const int pos_state = ++position & pos_state_mask;
    int next_price = cur_trial.price + price0( bm_match[cur_trial.state()][pos_state] );
    if( cur_trial.state.is_char() )
      next_price += literal_encoder.price_symbol( matchfinder[-1], cur_byte );
    else
      next_price += literal_encoder.price_matched( matchfinder[-1], cur_byte, match_byte );
    if( !matchfinder.move_pos() ) return 0;

    Trial & next_trial = trials[cur+1];

    next_trial.update( -1, cur, next_price );

    const int match_price = cur_trial.price + price1( bm_match[cur_trial.state()][pos_state] );
    const int rep_match_price = match_price + price1( bm_rep[cur_trial.state()] );

    if( match_byte == cur_byte && next_trial.dis != 0 )
      next_trial.update( 0, cur, rep_match_price +
                                 price_rep_len1( cur_trial.state, pos_state ) );

    const int len_limit = std::min( std::min( max_num_trials - 1 - cur,
                          matchfinder.available_bytes() ), match_len_limit );
    if( len_limit < min_match_len ) continue;

    for( int rep = 0; rep < num_rep_distances; ++rep )
      {
      const int dis = cur_trial.reps[rep] + 1;
      int len = 0;
      const uint8_t * const data = matchfinder.ptr_to_current_pos() - 1;
      while( len < len_limit && data[len] == data[len-dis] ) ++len;
      if( len >= min_match_len )
        {
        while( num_trials < cur + len ) trials[++num_trials].price = infinite_price;
        for( ; len >= min_match_len; --len )
          trials[cur+len].update( rep, cur, rep_match_price +
                                  price_rep( rep, len, cur_trial.state, pos_state ) );
        }
      }

    if( newlen <= len_limit &&
        ( newlen > min_match_len ||
          ( newlen == min_match_len && match_distances[newlen] < 256 ) ) )
      {
      const int normal_match_price = match_price + price0( bm_rep[cur_trial.state()] );
      while( num_trials < cur + newlen ) trials[++num_trials].price = infinite_price;

      for( int len = newlen; len >= min_match_len; --len )
        trials[cur+len].update( match_distances[len] + num_rep_distances, cur,
                                normal_match_price +
                                price_pair( match_distances[len], len, pos_state ) );
      }
    }
  }


     // End Of Stream mark => (dis == 0xFFFFFFFF, len == min_match_len)
void LZ_encoder::flush( const State & state )
  {
  const int pos_state = ( matchfinder.file_position() ) & pos_state_mask;
  range_encoder.encode_bit( bm_match[state()][pos_state], 1 );
  range_encoder.encode_bit( bm_rep[state()], 0 );
  const int len = min_match_len;
  len_encoder.encode( range_encoder, len - min_match_len, pos_state );
  const unsigned int dis = 0xFFFFFFFF;
  const int dis_slot = get_dis_slot( dis );
  range_encoder.encode_tree( bm_dis_slot[get_dis_state(len)], dis_slot, dis_slot_bits );
  const int direct_bits = ( dis_slot >> 1 ) - 1;
  const int direct_dis = dis & ( ( 1 << direct_bits ) - 1 );
  range_encoder.encode( direct_dis >> dis_align_bits, direct_bits - dis_align_bits );
  range_encoder.encode_tree_reversed( bm_align, direct_dis, dis_align_bits );
  range_encoder.flush_data();
  File_trailer trailer;
  trailer.file_crc( matchfinder.crc() );
  trailer.file_size( matchfinder.file_position() );
  for( unsigned int i = 0; i < sizeof trailer; ++i )
    range_encoder.put_byte( (( uint8_t *)&trailer)[i] );
  range_encoder.flush();
  }


LZ_encoder::LZ_encoder( const File_header & header, const int ides,
                        const int odes, const int len_limit )
  :
  match_len_limit( len_limit ),
  num_dis_slots( 2 * header.dictionary_bits ),
  longest_match_found( 0 ),
  matchfinder( header.dictionary_bits, match_len_limit, ides ),
  range_encoder( sizeof header, odes ),
  len_encoder( match_len_limit ),
  rep_match_len_encoder( match_len_limit ),
  literal_encoder()
  {
  if( header.dictionary_bits < min_dictionary_bits ||
      header.dictionary_bits > max_dictionary_bits ||
      match_len_limit < 5 || match_len_limit > max_match_len )
    internal_error( "invalid argument to encoder" );

  for( int slot = 0; slot < 4; ++slot ) dis_slots[slot] = slot;
  for( int c = 4, slot = 4; slot < 20; ++slot )
    {
    const int k = 1 << ( ( slot / 2 ) - 1 );
    for( int i = 0; i < k; ++i, ++c ) dis_slots[c] = slot;
    }
  fill_align_prices();
  }


bool LZ_encoder::encode()
  {
  int fill_counter = 0;
  int rep_distances[num_rep_distances];
  State state;
  uint8_t prev_byte = 0;
  for( int i = 0; i < num_rep_distances; ++i ) rep_distances[i] = 0;

  while( true )
    {
    if( matchfinder.finished() ) { flush( state ); return true; }
    if( fill_counter <= 0 ) { fill_distance_prices(); fill_counter = 512; }

    int ahead = best_pair_sequence( rep_distances, state );
    if( ahead <= 0 ) return false;
    fill_counter -= ahead;

    for( int i = 0; ahead > 0; )
      {
      const int pos_state = ( matchfinder.file_position() - ahead ) & pos_state_mask;
      int dis = trials[i].dis;
      const int len = trials[i].price;

      bool bit = ( dis < 0 && len == 1 );
      range_encoder.encode_bit( bm_match[state()][pos_state], !bit );
      if( bit )
        {
        const uint8_t cur_byte = matchfinder[-ahead];
        if( state.is_char() )
          literal_encoder.encode( range_encoder, prev_byte, cur_byte );
        else
          {
          const uint8_t match_byte = matchfinder[-rep_distances[0]-1-ahead];
          literal_encoder.encode_matched( range_encoder, prev_byte, match_byte, cur_byte );
          }
        state.set_char();
        prev_byte = cur_byte;
        }
      else
        {
        mtf_reps( dis, rep_distances );
        bit = ( dis < num_rep_distances );
        range_encoder.encode_bit( bm_rep[state()], bit );
        if( bit )
          {
          bit = ( dis == 0 );
          range_encoder.encode_bit( bm_rep0[state()], !bit );
          if( bit )
            range_encoder.encode_bit( bm_len[state()][pos_state], len > 1 );
          else
            {
            range_encoder.encode_bit( bm_rep1[state()], dis > 1 );
            if( dis > 1 )
              range_encoder.encode_bit( bm_rep2[state()], dis > 2 );
            }
          if( len == 1 ) state.set_short_rep();
          else
            {
            rep_match_len_encoder.encode( range_encoder, len - min_match_len, pos_state );
            state.set_rep();
            }
          }
        else
          {
          dis -= num_rep_distances;
          len_encoder.encode( range_encoder, len - min_match_len, pos_state );
          state.set_match();
          const int dis_slot = get_dis_slot( dis );
          range_encoder.encode_tree( bm_dis_slot[get_dis_state(len)], dis_slot, dis_slot_bits );

          if( dis_slot >= start_dis_model )
            {
            const int direct_bits = ( dis_slot >> 1 ) - 1;
            const int base = ( 2 | ( dis_slot & 1 ) ) << direct_bits;
            const int direct_dis = dis - base;

            if( dis_slot < end_dis_model )
              range_encoder.encode_tree_reversed( bm_dis + base - dis_slot,
                  direct_dis, direct_bits );
            else
              {
              range_encoder.encode( direct_dis >> dis_align_bits, direct_bits - dis_align_bits );
              range_encoder.encode_tree_reversed( bm_align, direct_dis, dis_align_bits );
              if( --align_price_count <= 0 ) fill_align_prices();
              }
            }
          }
        prev_byte = matchfinder[len-1-ahead];
        }
      ahead -= len; i += len;
      }
    }
  }
