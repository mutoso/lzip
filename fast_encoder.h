/*  Lzip - Data compressor based on the LZMA algorithm
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

class Fmatchfinder : public Matchfinder_base
  {
  enum { before_size = max_match_len + 1,
         dict_size = 65536,
         // bytes to keep in buffer after pos
         after_size = max_match_len,
         dict_factor = 16,
         len_limit = 16,
         num_prev_positions23 = 0,
         pos_array_factor = 1 };

  int key4;			// key made from latest 4 bytes

public:
  explicit Fmatchfinder( const int ifd )
    :
    Matchfinder_base( before_size, dict_size, after_size, dict_factor,
                      len_limit, num_prev_positions23, pos_array_factor, ifd ),
    key4( 0 )
    {}

  void reset() { Matchfinder_base::reset(); key4 = 0; }
  int longest_match_len( int * const distance );
  void longest_match_len( int n );
  };


class FLZ_encoder : public LZ_encoder_base
  {
  Fmatchfinder & fmatchfinder;

  void move_pos( int n )
    {
    if( --n >= 0 ) fmatchfinder.move_pos();
    fmatchfinder.longest_match_len( n );
    }

public:
  FLZ_encoder( Fmatchfinder & mf, const File_header & header, const int outfd )
    :
    LZ_encoder_base( header, mf.dictionary_size(), mf.match_len_limit(), outfd ),
    fmatchfinder( mf )
    {}

  bool encode_member( const unsigned long long member_size );
  };
