/*  Lzip - A data compressor based on the LZMA algorithm
    Copyright (C) 2008, 2009, 2010 Antonio Diaz Diaz.

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
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <stdint.h>

#include "lzip.h"
#include "decoder.h"


const CRC32 crc32;


bool Range_decoder::read_block()
  {
  if( !at_stream_end )
    {
    stream_pos = readblock( infd_, buffer, buffer_size );
    if( stream_pos != buffer_size && errno ) throw Error( "read error" );
    at_stream_end = ( stream_pos < buffer_size );
    partial_member_pos += pos;
    pos = 0;
    }
  return !finished();
  }


void LZ_decoder::flush_data()
  {
  const int size = pos - stream_pos;
  if( size > 0 )
    {
    crc32.update( crc_, buffer + stream_pos, size );
    if( outfd_ >= 0 &&
        writeblock( outfd_, buffer + stream_pos, size ) != size )
      throw Error( "write error" );
    if( pos >= buffer_size ) { partial_data_pos += pos; pos = 0; }
    stream_pos = pos;
    }
  }


bool LZ_decoder::verify_trailer( const Pretty_print & pp ) const
  {
  File_trailer trailer;
  const int trailer_size = File_trailer::size( member_version );
  const long long member_size = member_position() + trailer_size;
  bool error = false;

  for( int i = 0; i < trailer_size && !error; ++i )
    {
    if( !range_decoder.finished() )
      trailer.data[i] = range_decoder.get_byte();
    else
      {
      error = true;
      if( verbosity >= 0 )
        {
        pp();
        std::fprintf( stderr, "trailer truncated at trailer position %d;"
                              " some checks may fail.\n", i );
        }
      for( ; i < trailer_size; ++i ) trailer.data[i] = 0;
      }
    }
  if( member_version == 0 ) trailer.member_size( member_size );
  if( !range_decoder.code_is_zero() )
    {
    error = true;
    if( verbosity >= 0 )
      {
      pp();
      std::fprintf( stderr, "range_decoder final code is not zero.\n" );
      }
    }
  if( trailer.data_crc() != crc() )
    {
    error = true;
    if( verbosity >= 0 )
      {
      pp();
      std::fprintf( stderr, "crc mismatch; trailer says %08X, data crc is %08X.\n",
                    (unsigned int)trailer.data_crc(), (unsigned int)crc() );
      }
    }
  if( trailer.data_size() != data_position() )
    {
    error = true;
    if( verbosity >= 0 )
      {
      pp();
      std::fprintf( stderr, "data size mismatch; trailer says %lld, data size is %lld (0x%llX).\n",
                    trailer.data_size(), data_position(), data_position() );
      }
    }
  if( trailer.member_size() != member_size )
    {
    error = true;
    if( verbosity >= 0 )
      {
      pp();
      std::fprintf( stderr, "member size mismatch; trailer says %lld, member size is %lld (0x%llX).\n",
                    trailer.member_size(), member_size, member_size );
      }
    }
  if( !error && verbosity >= 3 )
    std::fprintf( stderr, "data crc %08X, data size %9lld, member size %8lld.  ",
                  (unsigned int)trailer.data_crc(), trailer.data_size(),
                  trailer.member_size() );
  return !error;
  }


    // Return value: 0 = OK, 1 = decoder error, 2 = unexpected EOF,
    //               3 = trailer error, 4 = unknown marker found.
int LZ_decoder::decode_member( const Pretty_print & pp )
  {
  unsigned int rep0 = 0;	// rep[0-3] latest four distances
  unsigned int rep1 = 0;	// used for efficient coding of
  unsigned int rep2 = 0;	// repeated distances
  unsigned int rep3 = 0;
  State state;
  range_decoder.load();

  while( true )
    {
    if( range_decoder.finished() ) { flush_data(); return 2; }
    const int pos_state = data_position() & pos_state_mask;
    if( range_decoder.decode_bit( bm_match[state()][pos_state] ) == 0 )
      {
      const uint8_t prev_byte = get_byte( 0 );
      if( state.is_char() )
        put_byte( literal_decoder.decode( range_decoder, prev_byte ) );
      else
        put_byte( literal_decoder.decode_matched( range_decoder, prev_byte,
                                                  get_byte( rep0 ) ) );
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
          unsigned int distance;
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
        unsigned int rep0_saved = rep0;
        len = min_match_len + len_decoder.decode( range_decoder, pos_state );
        const int dis_slot = range_decoder.decode_tree( bm_dis_slot[get_dis_state(len)], dis_slot_bits );
        if( dis_slot < start_dis_model ) rep0 = dis_slot;
        else
          {
          const int direct_bits = ( dis_slot >> 1 ) - 1;
          rep0 = ( 2 | ( dis_slot & 1 ) ) << direct_bits;
          if( dis_slot < end_dis_model )
            rep0 += range_decoder.decode_tree_reversed( bm_dis + rep0 - dis_slot, direct_bits );
          else
            {
            rep0 += range_decoder.decode( direct_bits - dis_align_bits ) << dis_align_bits;
            rep0 += range_decoder.decode_tree_reversed( bm_align, dis_align_bits );
            if( rep0 == 0xFFFFFFFFU )		// Marker found
              {
              rep0 = rep0_saved;
              range_decoder.normalize();
              flush_data();
              if( len == min_match_len )	// End Of Stream marker
                {
                if( verify_trailer( pp ) ) return 0; else return 3;
                }
              if( len == min_match_len + 1 )	// Sync Flush marker
                {
                range_decoder.load(); continue;
                }
              if( verbosity >= 0 )
                {
                pp();
                std::fprintf( stderr, "unsupported marker code `%d'.\n", len );
                }
              return 4;
              }
            if( rep0 >= (unsigned int)dictionary_size )
              { flush_data(); return 1; }
            }
          }
        rep3 = rep2; rep2 = rep1; rep1 = rep0_saved;
        state.set_match();
        }
      copy_block( rep0, len );
      }
    }
  }
