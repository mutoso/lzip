/*  Lzip - LZMA lossless data compressor
    Copyright (C) 2008, 2009, 2010, 2011, 2012, 2013 Antonio Diaz Diaz.

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


Dis_slots dis_slots;
Prob_prices prob_prices;


bool Matchfinder_base::read_block()
  {
  if( !at_stream_end && stream_pos < buffer_size )
    {
    const int size = buffer_size - stream_pos;
    const int rd = readblock( infd, buffer + stream_pos, size );
    stream_pos += rd;
    if( rd != size && errno ) throw Error( "Read error" );
    if( rd < size ) { at_stream_end = true; pos_limit = buffer_size; }
    }
  return pos < stream_pos;
  }


void Matchfinder_base::normalize_pos()
  {
  if( pos > stream_pos )
    internal_error( "pos > stream_pos in Matchfinder_base::normalize_pos" );
  if( !at_stream_end )
    {
    const int offset = pos - dictionary_size_ - before_size;
    const int size = stream_pos - offset;
    std::memmove( buffer, buffer + offset, size );
    partial_data_pos += offset;
    pos -= offset;
    stream_pos -= offset;
    for( int i = 0; i < num_prev_positions; ++i )
      if( prev_positions[i] >= 0 ) prev_positions[i] -= offset;
    for( int i = 0; i < pos_array_size; ++i )
      if( pos_array[i] >= 0 ) pos_array[i] -= offset;
    read_block();
    }
  }


Matchfinder_base::Matchfinder_base( const int before, const int dict_size,
                    const int after_size, const int dict_factor,
                    const int match_len_limit, const int num_prev_positions23,
                    const int pos_array_factor, const int ifd )
  :
  partial_data_pos( 0 ),
  before_size( before ),
  match_len_limit_( match_len_limit ),
  pos( 0 ),
  cyclic_pos( 0 ),
  stream_pos( 0 ),
  infd( ifd ),
  at_stream_end( false )
  {
  const int buffer_size_limit =
    ( dict_factor * dict_size ) + before_size + after_size;
  buffer_size = std::max( 65536, dict_size );
  buffer = (uint8_t *)std::malloc( buffer_size );
  if( !buffer ) throw std::bad_alloc();
  if( read_block() && !at_stream_end && buffer_size < buffer_size_limit )
    {
    buffer_size = buffer_size_limit;
    uint8_t * const tmp = (uint8_t *)std::realloc( buffer, buffer_size );
    if( !tmp ) { std::free( buffer ); throw std::bad_alloc(); }
    buffer = tmp;
    read_block();
    }
  if( at_stream_end && stream_pos < dict_size )
    dictionary_size_ = std::max( (int)min_dictionary_size, stream_pos );
  else
    dictionary_size_ = dict_size;
  pos_limit = buffer_size;
  if( !at_stream_end ) pos_limit -= after_size;
  unsigned size = 1 << std::max( 16, real_bits( dictionary_size_ - 1 ) - 2 );
  if( dictionary_size_ > 1 << 26 )		// 64 MiB
    size >>= 1;
  key4_mask = size - 1;
  size += num_prev_positions23;

  num_prev_positions = size;
  pos_array_size = pos_array_factor * ( dictionary_size_ + 1 );
  size += pos_array_size;
  prev_positions = new( std::nothrow ) int32_t[size];
  if( !prev_positions ) { std::free( buffer ); throw std::bad_alloc(); }
  pos_array = prev_positions + num_prev_positions;
  for( int i = 0; i < num_prev_positions; ++i ) prev_positions[i] = -1;
  }


void Matchfinder_base::reset()
  {
  const int size = stream_pos - pos;
  if( size > 0 ) std::memmove( buffer, buffer + pos, size );
  partial_data_pos = 0;
  stream_pos -= pos;
  pos = 0;
  cyclic_pos = 0;
  for( int i = 0; i < num_prev_positions; ++i ) prev_positions[i] = -1;
  read_block();
  }


int Matchfinder::get_match_pairs( struct Pair * pairs )
  {
  int len_limit = match_len_limit_;
  if( len_limit > available_bytes() )
    {
    len_limit = available_bytes();
    if( len_limit < 4 ) return 0;
    }

  int maxlen = min_match_len - 1;
  int num_pairs = 0;
  const int min_pos = ( pos > dictionary_size_ ) ? pos - dictionary_size_ : 0;
  const uint8_t * const data = buffer + pos;

  unsigned tmp = crc32[data[0]] ^ data[1];
  const int key2 = tmp & ( num_prev_positions2 - 1 );
  tmp ^= (uint32_t)data[2] << 8;
  const int key3 = num_prev_positions2 + ( tmp & ( num_prev_positions3 - 1 ) );
  const int key4 = num_prev_positions2 + num_prev_positions3 +
                   ( ( tmp ^ ( crc32[data[3]] << 5 ) ) & key4_mask );

  if( pairs )
    {
    int np2 = prev_positions[key2];
    int np3 = prev_positions[key3];
    if( np2 >= min_pos && buffer[np2] == data[0] )
      {
      pairs[0].dis = pos - np2 - 1;
      pairs[0].len = maxlen = 2;
      num_pairs = 1;
      }
    if( np2 != np3 && np3 >= min_pos && buffer[np3] == data[0] )
      {
      maxlen = 3;
      pairs[num_pairs].dis = pos - np3 - 1;
      ++num_pairs;
      np2 = np3;
      }
    if( num_pairs > 0 )
      {
      const int delta = pos - np2;
      while( maxlen < len_limit && data[maxlen-delta] == data[maxlen] )
        ++maxlen;
      pairs[num_pairs-1].len = maxlen;
      if( maxlen >= len_limit ) pairs = 0;
      }
    if( maxlen < 3 ) maxlen = 3;
    }

  prev_positions[key2] = pos;
  prev_positions[key3] = pos;
  int newpos = prev_positions[key4];
  prev_positions[key4] = pos;

  int32_t * ptr0 = pos_array + ( cyclic_pos << 1 );
  int32_t * ptr1 = ptr0 + 1;
  int len = 0, len0 = 0, len1 = 0;

  for( int count = cycles; ; )
    {
    if( newpos < min_pos || --count < 0 ) { *ptr0 = *ptr1 = -1; break; }

    const int delta = pos - newpos;
    int32_t * const newptr = pos_array +
      ( ( cyclic_pos - delta +
          ( ( cyclic_pos >= delta ) ? 0 : dictionary_size_ + 1 ) ) << 1 );
    if( data[len-delta] == data[len] )
      {
      while( ++len < len_limit && data[len-delta] == data[len] ) {}
      if( pairs && maxlen < len )
        {
        pairs[num_pairs].dis = delta - 1;
        pairs[num_pairs].len = maxlen = len;
        ++num_pairs;
        }
      if( len >= len_limit )
        {
        *ptr0 = newptr[0];
        *ptr1 = newptr[1];
        break;
        }
      }
    if( data[len-delta] < data[len] )
      {
      *ptr0 = newpos;
      ptr0 = newptr + 1;
      newpos = *ptr0;
      len0 = len; if( len1 < len ) len = len1;
      }
    else
      {
      *ptr1 = newpos;
      ptr1 = newptr;
      newpos = *ptr1;
      len1 = len; if( len0 < len ) len = len0;
      }
    }
  return num_pairs;
  }


void Range_encoder::flush_data()
  {
  if( pos > 0 )
    {
    if( outfd >= 0 && writeblock( outfd, buffer, pos ) != pos )
      throw Error( "Write error" );
    partial_member_pos += pos;
    pos = 0;
    if( verbosity >= 2 ) show_progress();
    }
  }


void Len_encoder::encode( Range_encoder & renc, int symbol,
                          const int pos_state )
  {
  symbol -= min_match_len;
  if( symbol < len_low_symbols )
    {
    renc.encode_bit( choice1, 0 );
    renc.encode_tree( bm_low[pos_state], symbol, len_low_bits );
    }
  else
    {
    renc.encode_bit( choice1, 1 );
    if( symbol < len_low_symbols + len_mid_symbols )
      {
      renc.encode_bit( choice2, 0 );
      renc.encode_tree( bm_mid[pos_state], symbol - len_low_symbols,
                        len_mid_bits );
      }
    else
      {
      renc.encode_bit( choice2, 1 );
      renc.encode_tree( bm_high, symbol - len_low_symbols - len_mid_symbols,
                        len_high_bits );
      }
    }
  if( --counters[pos_state] <= 0 ) update_prices( pos_state );
  }


     // End Of Stream mark => (dis == 0xFFFFFFFFU, len == min_match_len)
void LZ_encoder_base::full_flush( const unsigned long long data_position,
                                  const State state )
  {
  const int pos_state = data_position & pos_state_mask;
  renc.encode_bit( bm_match[state()][pos_state], 1 );
  renc.encode_bit( bm_rep[state()], 0 );
  encode_pair( 0xFFFFFFFFU, min_match_len, pos_state );
  renc.flush();
  File_trailer trailer;
  trailer.data_crc( crc() );
  trailer.data_size( data_position );
  trailer.member_size( renc.member_position() + File_trailer::size() );
  for( int i = 0; i < File_trailer::size(); ++i )
    renc.put_byte( trailer.data[i] );
  renc.flush_data();
  }


void LZ_encoder::fill_align_prices()
  {
  for( int i = 0; i < dis_align_size; ++i )
    align_prices[i] = price_symbol_reversed( bm_align, i, dis_align_bits );
  align_price_count = dis_align_size;
  }


void LZ_encoder::fill_distance_prices()
  {
  for( int dis = start_dis_model; dis < modeled_distances; ++dis )
    {
    const int dis_slot = dis_slots[dis];
    const int direct_bits = ( dis_slot >> 1 ) - 1;
    const int base = ( 2 | ( dis_slot & 1 ) ) << direct_bits;
    const int price = price_symbol_reversed( bm_dis + base - dis_slot - 1,
                                             dis - base, direct_bits );
    for( int len_state = 0; len_state < len_states; ++len_state )
      dis_prices[len_state][dis] = price;
    }

  for( int len_state = 0; len_state < len_states; ++len_state )
    {
    int * const dsp = dis_slot_prices[len_state];
    const Bit_model * const bmds = bm_dis_slot[len_state];
    int slot = 0;
    for( ; slot < end_dis_model && slot < num_dis_slots; ++slot )
      dsp[slot] = price_symbol( bmds, slot, dis_slot_bits );
    for( ; slot < num_dis_slots; ++slot )
      dsp[slot] = price_symbol( bmds, slot, dis_slot_bits ) +
                  (((( slot >> 1 ) - 1 ) - dis_align_bits ) << price_shift_bits );

    int * const dp = dis_prices[len_state];
    int dis = 0;
    for( ; dis < start_dis_model; ++dis )
      dp[dis] = dsp[dis];
    for( ; dis < modeled_distances; ++dis )
      dp[dis] += dsp[dis_slots[dis]];
    }
  }


/* Return value == number of bytes advanced (ahead).
   trials[0]..trials[ahead-1] contain the steps to encode.
   ( trials[0].dis == -1 && trials[0].price == 1 ) means literal.
*/
int LZ_encoder::sequence_optimizer( const int reps[num_rep_distances],
                                    const State state )
  {
  int num_pairs, num_trials;

  if( pending_num_pairs > 0 )			// from previous call
    {
    num_pairs = pending_num_pairs;
    pending_num_pairs = 0;
    }
  else
    num_pairs = read_match_distances();
  const int main_len = ( num_pairs > 0 ) ? pairs[num_pairs-1].len : 0;

  int replens[num_rep_distances];
  int rep_index = 0;
  for( int i = 0; i < num_rep_distances; ++i )
    {
    replens[i] = matchfinder.true_match_len( 0, reps[i] + 1, max_match_len );
    if( replens[i] > replens[rep_index] ) rep_index = i;
    }
  if( replens[rep_index] >= matchfinder.match_len_limit() )
    {
    trials[0].dis = rep_index;
    trials[0].price = replens[rep_index];
    move_pos( replens[rep_index] );
    return replens[rep_index];
    }

  if( main_len >= matchfinder.match_len_limit() )
    {
    trials[0].dis = pairs[num_pairs-1].dis + num_rep_distances;
    trials[0].price = main_len;
    move_pos( main_len );
    return main_len;
    }

  const int pos_state = matchfinder.data_position() & pos_state_mask;
  const uint8_t prev_byte = matchfinder[-1];
  const uint8_t cur_byte = matchfinder[0];
  const uint8_t match_byte = matchfinder[-reps[0]-1];

  trials[0].state = state;
  trials[1].dis = -1;
  trials[1].price = price0( bm_match[state()][pos_state] );
  if( state.is_char() )
    trials[1].price += price_literal( prev_byte, cur_byte );
  else
    trials[1].price += price_matched( prev_byte, cur_byte, match_byte );

  const int match_price = price1( bm_match[state()][pos_state] );
  const int rep_match_price = match_price + price1( bm_rep[state()] );

  if( match_byte == cur_byte )
    trials[1].update( rep_match_price + price_rep_len1( state, pos_state ), 0, 0 );

  num_trials = std::max( main_len, replens[rep_index] );

  if( num_trials < min_match_len )
    {
    trials[0].dis = trials[1].dis;
    trials[0].price = 1;
    matchfinder.move_pos();
    return 1;
    }

  for( int i = 0; i < num_rep_distances; ++i )
    trials[0].reps[i] = reps[i];
  trials[1].prev_index = 0;
  trials[1].prev_index2 = single_step_trial;

  for( int len = min_match_len; len <= num_trials; ++len )
    trials[len].price = infinite_price;

  for( int rep = 0; rep < num_rep_distances; ++rep )
    {
    if( replens[rep] < min_match_len ) continue;
    const int price = rep_match_price + price_rep( rep, state, pos_state );
    for( int len = min_match_len; len <= replens[rep]; ++len )
      trials[len].update( price + rep_len_encoder.price( len, pos_state ),
                          rep, 0 );
    }

  if( main_len > replens[0] )
    {
    const int normal_match_price = match_price + price0( bm_rep[state()] );
    int i = 0, len = std::max( replens[0] + 1, (int)min_match_len );
    while( len > pairs[i].len ) ++i;
    while( true )
      {
      const int dis = pairs[i].dis;
      trials[len].update( normal_match_price + price_pair( dis, len, pos_state ),
                          dis + num_rep_distances, 0 );
      if( ++len > pairs[i].len && ++i >= num_pairs ) break;
      }
    }

  int cur = 0;
  matchfinder.move_pos();

  while( true )				// price optimization loop
    {
    if( ++cur >= num_trials )		// no more initialized trials
      {
      backward( cur );
      return cur;
      }

    const int num_pairs = read_match_distances();
    const int newlen = ( num_pairs > 0 ) ? pairs[num_pairs-1].len : 0;
    if( newlen >= matchfinder.match_len_limit() )
      {
      pending_num_pairs = num_pairs;
      backward( cur );
      return cur;
      }

    // give final values to current trial
    Trial & cur_trial = trials[cur];
    int prev_index = cur_trial.prev_index;
    const int prev_index2 = cur_trial.prev_index2;
    State cur_state;

    if( prev_index2 != single_step_trial )
      {
      --prev_index;
      if( prev_index2 >= 0 )
        {
        cur_state = trials[prev_index2].state;
        if( cur_trial.dis2 < num_rep_distances )
          cur_state.set_rep();
        else
          cur_state.set_match();
        }
      else
        cur_state = trials[prev_index].state;
      cur_state.set_char();
      }
    else
      cur_state = trials[prev_index].state;

    if( prev_index == cur - 1 )
      {
      if( cur_trial.dis == 0 ) cur_state.set_short_rep();
      else cur_state.set_char();
      for( int i = 0; i < num_rep_distances; ++i )
        cur_trial.reps[i] = trials[prev_index].reps[i];
      }
    else
      {
      int dis;
      if( prev_index2 >= 0 )
        {
        dis = cur_trial.dis2;
        prev_index = prev_index2;
        cur_state.set_rep();
        }
      else
        {
        dis = cur_trial.dis;
        if( dis < num_rep_distances ) cur_state.set_rep();
        else cur_state.set_match();
        }
      for( int i = 0; i < num_rep_distances; ++i )
        cur_trial.reps[i] = trials[prev_index].reps[i];
      mtf_reps( dis, cur_trial.reps );
      }
    cur_trial.state = cur_state;

    const int pos_state = matchfinder.data_position() & pos_state_mask;
    const uint8_t prev_byte = matchfinder[-1];
    const uint8_t cur_byte = matchfinder[0];
    const uint8_t match_byte = matchfinder[-cur_trial.reps[0]-1];

    int next_price = cur_trial.price +
                     price0( bm_match[cur_state()][pos_state] );
    if( cur_state.is_char() )
      next_price += price_literal( prev_byte, cur_byte );
    else
      next_price += price_matched( prev_byte, cur_byte, match_byte );
    matchfinder.move_pos();

    // try last updates to next trial
    Trial & next_trial = trials[cur+1];

    next_trial.update( next_price, -1, cur );

    const int match_price = cur_trial.price + price1( bm_match[cur_state()][pos_state] );
    const int rep_match_price = match_price + price1( bm_rep[cur_state()] );

    if( match_byte == cur_byte && next_trial.dis != 0 )
      {
      const int price = rep_match_price +
                        price_rep_len1( cur_state, pos_state );
      if( price <= next_trial.price )
        {
        next_trial.price = price;
        next_trial.dis = 0;
        next_trial.prev_index = cur;
        next_trial.prev_index2 = single_step_trial;
        }
      }

    const int available_bytes = std::min( matchfinder.available_bytes() + 1,
                                          max_num_trials - 1 - cur );
    if( available_bytes < min_match_len ) continue;

    const int len_limit = std::min( matchfinder.match_len_limit(),
                                    available_bytes );

    // try literal + rep0
    if( match_byte != cur_byte && next_trial.prev_index != cur )
      {
      const uint8_t * const data = matchfinder.ptr_to_current_pos() - 1;
      const int dis = cur_trial.reps[0] + 1;
      const int limit = std::min( matchfinder.match_len_limit() + 1,
                                  available_bytes );
      int len = 1;
      while( len < limit && data[len-dis] == data[len] ) ++len;
      if( --len >= min_match_len )
        {
        const int pos_state2 = ( pos_state + 1 ) & pos_state_mask;
        State state2 = cur_state; state2.set_char();
        const int price = next_price +
                  price1( bm_match[state2()][pos_state2] ) +
                  price1( bm_rep[state2()] ) +
                  price_rep0_len( len, state2, pos_state2 );
        while( num_trials < cur + 1 + len )
          trials[++num_trials].price = infinite_price;
        trials[cur+1+len].update2( price, 0, cur + 1 );
        }
      }

    int start_len = min_match_len;

    // try rep distances
    for( int rep = 0; rep < num_rep_distances; ++rep )
      {
      const uint8_t * const data = matchfinder.ptr_to_current_pos() - 1;
      int len;
      const int dis = cur_trial.reps[rep] + 1;

      if( data[-dis] != data[0] || data[1-dis] != data[1] ) continue;
      for( len = min_match_len; len < len_limit; ++len )
        if( data[len-dis] != data[len] ) break;
      while( num_trials < cur + len )
        trials[++num_trials].price = infinite_price;
      int price = rep_match_price + price_rep( rep, cur_state, pos_state );
      for( int i = min_match_len; i <= len; ++i )
        trials[cur+i].update( price + rep_len_encoder.price( i, pos_state ),
                              rep, cur );

      if( rep == 0 ) start_len = len + 1;	// discard shorter matches

      // try rep + literal + rep0
      int len2 = len + 1;
      const int limit = std::min( matchfinder.match_len_limit() + len2,
                                  available_bytes );
      while( len2 < limit && data[len2-dis] == data[len2] ) ++len2;
      len2 -= len + 1;
      if( len2 < min_match_len ) continue;

      int pos_state2 = ( pos_state + len ) & pos_state_mask;
      State state2 = cur_state; state2.set_rep();
      price += rep_len_encoder.price( len, pos_state ) +
               price0( bm_match[state2()][pos_state2] ) +
               price_matched( data[len-1], data[len], data[len-dis] );
      pos_state2 = ( pos_state2 + 1 ) & pos_state_mask;
      state2.set_char();
      price += price1( bm_match[state2()][pos_state2] ) +
               price1( bm_rep[state2()] ) +
               price_rep0_len( len2, state2, pos_state2 );
      while( num_trials < cur + len + 1 + len2 )
        trials[++num_trials].price = infinite_price;
      trials[cur+len+1+len2].update3( price, 0, cur + len + 1, rep, cur );
      }

    // try matches
    if( newlen >= start_len && newlen <= len_limit )
      {
      const int normal_match_price = match_price +
                                     price0( bm_rep[cur_state()] );

      while( num_trials < cur + newlen )
        trials[++num_trials].price = infinite_price;

      int i = 0;
      while( start_len > pairs[i].len ) ++i;
      int dis = pairs[i].dis;
      for( int len = start_len; ; ++len )
        {
        int price = normal_match_price + price_pair( dis, len, pos_state );

        trials[cur+len].update( price, dis + num_rep_distances, cur );

        // try match + literal + rep0
        if( len == pairs[i].len )
          {
          const uint8_t * const data = matchfinder.ptr_to_current_pos() - 1;
          const int dis2 = dis + 1;
          int len2 = len + 1;
          const int limit = std::min( matchfinder.match_len_limit() + len2,
                                      available_bytes );
          while( len2 < limit && data[len2-dis2] == data[len2] ) ++len2;
          len2 -= len + 1;
          if( len2 >= min_match_len )
            {
            int pos_state2 = ( pos_state + len ) & pos_state_mask;
            State state2 = cur_state; state2.set_match();
            price += price0( bm_match[state2()][pos_state2] ) +
                     price_matched( data[len-1], data[len], data[len-dis2] );
            pos_state2 = ( pos_state2 + 1 ) & pos_state_mask;
            state2.set_char();
            price += price1( bm_match[state2()][pos_state2] ) +
                     price1( bm_rep[state2()] ) +
                     price_rep0_len( len2, state2, pos_state2 );

            while( num_trials < cur + len + 1 + len2 )
              trials[++num_trials].price = infinite_price;
            trials[cur+len+1+len2].update3( price, 0, cur + len + 1,
                                            dis + num_rep_distances, cur );
            }
          if( ++i >= num_pairs ) break;
          dis = pairs[i].dis;
          }
        }
      }
    }
  }


bool LZ_encoder::encode_member( const unsigned long long member_size )
  {
  const unsigned long long member_size_limit =
    member_size - File_trailer::size() - max_marker_size;
  const int fill_count = ( matchfinder.match_len_limit() > 12 ) ? 128 : 512;
  int fill_counter = 0;
  int rep_distances[num_rep_distances];
  State state;
  for( int i = 0; i < num_rep_distances; ++i ) rep_distances[i] = 0;

  if( matchfinder.data_position() != 0 ||
      renc.member_position() != File_header::size )
    return false;			// can be called only once

  if( !matchfinder.finished() )		// encode first byte
    {
    const uint8_t prev_byte = 0;
    const uint8_t cur_byte = matchfinder[0];
    renc.encode_bit( bm_match[state()][0], 0 );
    encode_literal( prev_byte, cur_byte );
    crc32.update_byte( crc_, cur_byte );
    matchfinder.get_match_pairs();
    matchfinder.move_pos();
    }

  while( !matchfinder.finished() )
    {
    if( pending_num_pairs == 0 )
      {
      if( fill_counter <= 0 )
        { fill_distance_prices(); fill_counter = fill_count; }
      if( align_price_count <= 0 ) fill_align_prices();
      }

    int ahead = sequence_optimizer( rep_distances, state );
    if( ahead <= 0 ) return false;		// can't happen

    for( int i = 0; ; )
      {
      const int pos_state =
        ( matchfinder.data_position() - ahead ) & pos_state_mask;
      const int dis = trials[i].dis;
      const int len = trials[i].price;

      bool bit = ( dis < 0 && len == 1 );
      renc.encode_bit( bm_match[state()][pos_state], !bit );
      if( bit )					// literal byte
        {
        const uint8_t prev_byte = matchfinder[-ahead-1];
        const uint8_t cur_byte = matchfinder[-ahead];
        crc32.update_byte( crc_, cur_byte );
        if( state.is_char() )
          encode_literal( prev_byte, cur_byte );
        else
          {
          const uint8_t match_byte = matchfinder[-ahead-rep_distances[0]-1];
          encode_matched( prev_byte, cur_byte, match_byte );
          }
        state.set_char();
        }
      else					// match or repeated match
        {
        crc32.update_buf( crc_, matchfinder.ptr_to_current_pos() - ahead, len );
        mtf_reps( dis, rep_distances );
        bit = ( dis < num_rep_distances );
        renc.encode_bit( bm_rep[state()], bit );
        if( bit )
          {
          bit = ( dis == 0 );
          renc.encode_bit( bm_rep0[state()], !bit );
          if( bit )
            renc.encode_bit( bm_len[state()][pos_state], len > 1 );
          else
            {
            renc.encode_bit( bm_rep1[state()], dis > 1 );
            if( dis > 1 )
              renc.encode_bit( bm_rep2[state()], dis > 2 );
            }
          if( len == 1 ) state.set_short_rep();
          else
            {
            rep_len_encoder.encode( renc, len, pos_state );
            state.set_rep();
            }
          }
        else
          {
          encode_pair( dis - num_rep_distances, len, pos_state );
          if( get_slot( dis - num_rep_distances ) >= end_dis_model )
            --align_price_count;
          --fill_counter;
          state.set_match();
          }
        }
      ahead -= len; i += len;
      if( renc.member_position() >= member_size_limit )
        {
        if( !matchfinder.dec_pos( ahead ) ) return false;
        full_flush( matchfinder.data_position(), state );
        return true;
        }
      if( ahead <= 0 ) break;
      }
    }
  full_flush( matchfinder.data_position(), state );
  return true;
  }
