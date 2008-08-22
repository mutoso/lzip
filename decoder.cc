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

#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <stdint.h>

#include "lzip.h"
#include "decoder.h"


const Update_crc update_crc;


bool Input_buffer::read_block()
  {
  if( _finished ) return false;
  size = readblock( _ides, (char *)buffer, buffer_size );
  if( size != buffer_size && errno ) throw Error( "read error" );
  pos = 0;
  _finished = ( size == 0 );
  return !_finished;
  }


void LZ_decoder::flush()
  {
  if( stream_pos == 0 )
    {
    const int wr = writeblock( _odes, (char *)buffer, pos );
    if( wr != pos ) throw Error( "write error" );
    if( pos >= buffer_size ) { partial_file_pos += pos; pos = 0; }
    stream_pos = pos;
    }
  }


long long LZ_decoder::verify_trailer( const Pretty_print & pp )
  {
  flush();
  const long long file_size = file_position();
  const uint32_t file_crc = crc();
  bool error = false;
  File_trailer trailer;
  for( unsigned int i = 0; i < sizeof trailer; ++i )
    ((uint8_t *)&trailer)[i] = range_decoder.read_byte();
  if( trailer.file_crc() != file_crc )
    {
    if( verbosity >= 0 ) pp( "bad crc for uncompressed file" );
    error = true;
    }
  if( trailer.file_size() != file_size )
    {
    if( verbosity >= 0 )
      {
      if( trailer.file_size() >= 0 )
        { pp();
          std::fprintf( stderr, "bad uncompressed size; expected %lld, obtained %lld.\n",
                        trailer.file_size(), file_size ); }
      else pp( "bad file trailer" );
      }
    error = true;
    }
  if( error ) return range_decoder.file_position() - sizeof trailer;
  return -1;
  }


     // Returns position of failure or -1 if no error.
long long LZ_decoder::decode( const Pretty_print & pp )
  {
  int rep0 = 0;
  int rep1 = 0;
  int rep2 = 0;
  int rep3 = 0;
  State state;
  uint8_t prev_byte = 0;

  while( true )
    {
    if( range_decoder.finished() ) return range_decoder.file_position();
    const int pos_state = file_position() & pos_state_mask;
    if( range_decoder.decode_bit( bm_match[state()][pos_state] ) == 0 )
      {
      if( state.is_char() )
        prev_byte = literal_decoder.decode( range_decoder, prev_byte );
      else
        prev_byte = literal_decoder.decode_matched( range_decoder, prev_byte,
                                                    get_byte( rep0 ) );
      put_byte( prev_byte );
      state.set_char();
      }
    else
      {
      int len;
      if( range_decoder.decode_bit( bm_rep[state()] ) == 1 )
        {
        len = 0;
        if( range_decoder.decode_bit( bm_rep0[state()] ) == 0 )
          {
          if( range_decoder.decode_bit( bm_len[state()][pos_state] ) == 0 )
            { len = 1; state.set_short_rep(); }
          }
        else
          {
          int distance;
          if( range_decoder.decode_bit( bm_rep1[state()] ) == 0 )
            distance = rep1;
          else
            {
            if( range_decoder.decode_bit( bm_rep2[state()] ) == 0 )
              distance = rep2;
            else { distance = rep3; rep3 = rep2; }
            rep2 = rep1;
            }
          rep1 = rep0;
          rep0 = distance;
          }
        if( len == 0 )
          {
          len = min_match_len + rep_match_len_decoder.decode( range_decoder, pos_state );
          state.set_rep();
          }
        }
      else
        {
        rep3 = rep2; rep2 = rep1; rep1 = rep0;
        len = min_match_len + len_decoder.decode( range_decoder, pos_state );
        state.set_match();
        const int dis_slot = range_decoder.decode_tree( bm_dis_slot[get_dis_state(len)], dis_slot_bits );
        if( dis_slot < start_dis_model ) rep0 = dis_slot;
        else
          {
          const int direct_bits = ( dis_slot >> 1 ) - 1;
          unsigned int dis = ( 2 | ( dis_slot & 1 ) ) << direct_bits;
          if( dis_slot < end_dis_model )
            dis += range_decoder.decode_tree_reversed( bm_dis + dis - dis_slot, direct_bits );
          else
            {
            dis += range_decoder.decode( direct_bits - dis_align_bits ) << dis_align_bits;
            dis += range_decoder.decode_tree_reversed( bm_align, dis_align_bits );
            if( dis == 0xFFFFFFFF )	// End Of Stream mark
              {
              if( len == min_match_len ) return verify_trailer( pp );
              return range_decoder.file_position();
              }
            }
          rep0 = dis;
          }
        }
      if( !copy_block( rep0, len ) ) return range_decoder.file_position();
      prev_byte = get_byte( 0 );
      }
    }
  }
