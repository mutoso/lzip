/*  Lzip - LZMA lossless data compressor
    Copyright (C) 2008-2014 Antonio Diaz Diaz.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

class Fmatchfinder : public Matchfinder_base
  {
  enum { before_size = max_match_len + 1,
         dict_size = 65536,
         // bytes to keep in buffer after pos
         after_size = max_match_len,
         dict_factor = 16,
         num_prev_positions23 = 0,
         pos_array_factor = 1 };

  int key4;			// key made from latest 4 bytes

  void reset_key4()
    {
    key4 = 0;
    for( int i = 0; i < 3 && i < available_bytes(); ++i )
      key4 = ( key4 << 4 ) ^ buffer[i];
    }

public:
  explicit Fmatchfinder( const int ifd )
    :
    Matchfinder_base( before_size, dict_size, after_size, dict_factor,
                      num_prev_positions23, pos_array_factor, ifd )
    { reset_key4(); }

  enum { len_limit = 16 };
  void reset() { Matchfinder_base::reset(); reset_key4(); }
  int longest_match_len( int * const distance );

  void update_and_move( int n )
    {
    while( --n >= 0 )
      {
      if( available_bytes() >= 4 )
        {
        key4 = ( ( key4 << 4 ) ^ buffer[pos+3] ) & key4_mask;
        const int newpos = prev_positions[key4];
        prev_positions[key4] = pos + 1;
        pos_array[cyclic_pos] = newpos;
        }
      move_pos();
      }
    }
  };


class FLZ_encoder : public LZ_encoder_base
  {
  Fmatchfinder & fmatchfinder;

public:
  FLZ_encoder( Fmatchfinder & mf, const File_header & header, const int outfd )
    :
    LZ_encoder_base( header, outfd ),
    fmatchfinder( mf )
    {}

  bool encode_member( const unsigned long long member_size );
  };
