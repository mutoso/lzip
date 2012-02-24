/*  Lzip - Data compressor based on the LZMA algorithm
    Copyright (C) 2008, 2009, 2010, 2011, 2012 Antonio Diaz Diaz.

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
#include "fast_encoder.h"


int Fmatchfinder::longest_match_len( int * const distance )
  {
  int len_limit = match_len_limit_;
  if( len_limit > available_bytes() )
    {
    len_limit = available_bytes();
    if( len_limit < 4 ) return 0;
    }

  key4 = ( ( key4 << 4 ) ^ buffer[pos+3] ) & ( num_prev_pos - 1 );

  int newpos = prev_positions[key4];
  prev_positions[key4] = pos;

  int32_t * ptr0 = pos_array + cyclic_pos;
  int maxlen = 0;

  for( int count = 4; ; )
    {
    if( newpos < (pos - dictionary_size_ + 1) || newpos < 0 || --count < 0 )
      { *ptr0 = -1; break; }
    int len = 0;
    if( buffer[maxlen+newpos] == buffer[maxlen+pos] )
      while( len < len_limit && buffer[len+newpos] == buffer[len+pos] ) ++len;

    const int delta = pos - newpos;
    if( maxlen < len ) { maxlen = len; *distance = delta - 1; }

    int32_t * const newptr = pos_array +
      ( cyclic_pos - delta +
          ( ( cyclic_pos >= delta ) ? 0 : dictionary_size_ ) );

    if( len < len_limit )
      {
      *ptr0 = newpos;
      ptr0 = newptr;
      newpos = *ptr0;
      }
    else
      {
      *ptr0 = *newptr;
      break;
      }
    }
  if( maxlen == match_len_limit_ && maxlen < max_match_len )
    maxlen += true_match_len( maxlen, *distance + 1, max_match_len - maxlen );
  return maxlen;
  }


void Fmatchfinder::longest_match_len()
  {
  int len_limit = match_len_limit_;
  if( len_limit > available_bytes() )
    {
    len_limit = available_bytes();
    if( len_limit < 4 ) return;
    }

  key4 = ( ( key4 << 4 ) ^ buffer[pos+3] ) & ( num_prev_pos - 1 );

  const int newpos = prev_positions[key4];
  prev_positions[key4] = pos;

  int32_t * const ptr0 = pos_array + cyclic_pos;

  if( newpos < (pos - dictionary_size_ + 1) || newpos < 0 )
    *ptr0 = -1;
  else if( buffer[len_limit-1+newpos] != buffer[len_limit-1+pos] ||
           std::memcmp( buffer + newpos, buffer + pos, len_limit - 1 ) )
    *ptr0 = newpos;
  else
    {
    int idx = cyclic_pos - pos + newpos;
    if( idx < 0 ) idx += dictionary_size_;
    *ptr0 = pos_array[idx];
    }
  }


void FLZ_encoder::sequence_optimizer( int reps[num_rep_distances],
                                      State & state )
  {
  int match_distance;
  const int main_len = fmatchfinder.longest_match_len( &match_distance );
  const int pos_state = fmatchfinder.data_position() & pos_state_mask;
  int dis = 0;
  int len = 0;

  for( int i = 0; i < num_rep_distances; ++i )
    {
    const int tlen =
      fmatchfinder.true_match_len( 0, reps[i] + 1, max_match_len );
    if( tlen > len ) { len = tlen; dis = i; }
    }
  if( len > min_match_len && len + 4 > main_len )
    {
    crc32.update( crc_, fmatchfinder.ptr_to_current_pos(), len );
    mtf_reps( dis, reps );
    range_encoder.encode_bit( bm_match[state()][pos_state], 1 );
    range_encoder.encode_bit( bm_rep[state()], 1 );
    const bool bit = ( dis == 0 );
    range_encoder.encode_bit( bm_rep0[state()], !bit );
    if( bit )
      range_encoder.encode_bit( bm_len[state()][pos_state], 1 );
    else
      {
      range_encoder.encode_bit( bm_rep1[state()], dis > 1 );
      if( dis > 1 )
        range_encoder.encode_bit( bm_rep2[state()], dis > 2 );
      }
    rep_match_len_encoder.encode( range_encoder, len, pos_state );
    state.set_rep();
    move_pos( len );
    return;
    }

  if( main_len > min_match_len ||
      ( main_len == min_match_len && match_distance < modeled_distances ) )
    {
    crc32.update( crc_, fmatchfinder.ptr_to_current_pos(), main_len );
    dis = match_distance;
    mtf_reps( dis + num_rep_distances, reps );
    range_encoder.encode_bit( bm_match[state()][pos_state], 1 );
    range_encoder.encode_bit( bm_rep[state()], 0 );
    encode_pair( dis, main_len, pos_state );
    state.set_match();
    move_pos( main_len );
    return;
    }

  const uint8_t prev_byte = fmatchfinder[-1];
  const uint8_t cur_byte = fmatchfinder[0];
  const uint8_t match_byte = fmatchfinder[-reps[0]-1];
  crc32.update( crc_, cur_byte );
  fmatchfinder.move_pos();

  if( match_byte == cur_byte )
    {
    int price = price0( bm_match[state()][pos_state] );
    if( state.is_char() )
      price += literal_encoder.price_symbol( prev_byte, cur_byte );
    else
      price += literal_encoder.price_matched( prev_byte, cur_byte, match_byte );
    const int short_rep_price = price1( bm_match[state()][pos_state] ) +
                                price1( bm_rep[state()] ) +
                                price0( bm_rep0[state()] ) +
                                price0( bm_len[state()][pos_state] );
    if( short_rep_price < price )
      {
      range_encoder.encode_bit( bm_match[state()][pos_state], 1 );
      range_encoder.encode_bit( bm_rep[state()], 1 );
      range_encoder.encode_bit( bm_rep0[state()], 0 );
      range_encoder.encode_bit( bm_len[state()][pos_state], 0 );
      state.set_short_rep();
      return;
      }
    }

  // literal byte
  range_encoder.encode_bit( bm_match[state()][pos_state], 0 );
  if( state.is_char() )
    literal_encoder.encode( range_encoder, prev_byte, cur_byte );
  else
    literal_encoder.encode_matched( range_encoder,
                                    prev_byte, cur_byte, match_byte );
  state.set_char();
  }


bool FLZ_encoder::encode_member( const long long member_size )
  {
  const long long member_size_limit =
    member_size - File_trailer::size() - max_marker_size;
  int rep_distances[num_rep_distances];
  State state;
  for( int i = 0; i < num_rep_distances; ++i ) rep_distances[i] = 0;

  if( fmatchfinder.data_position() != 0 ||
      range_encoder.member_position() != File_header::size )
    return false;			// can be called only once

  if( !fmatchfinder.finished() )	// encode first byte
    {
    const uint8_t prev_byte = 0;
    const uint8_t cur_byte = fmatchfinder[0];
    range_encoder.encode_bit( bm_match[state()][0], 0 );
    literal_encoder.encode( range_encoder, prev_byte, cur_byte );
    crc32.update( crc_, cur_byte );
    fmatchfinder.longest_match_len();
    fmatchfinder.move_pos();
    }

  while( !fmatchfinder.finished() &&
         range_encoder.member_position() < member_size_limit )
    sequence_optimizer( rep_distances, state );

  full_flush( fmatchfinder.data_position(), state );
  return true;
  }
