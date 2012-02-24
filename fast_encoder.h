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

class Fmatchfinder : public Matchfinder_base
  {
  enum { before = max_match_len + 1,
         dict_size = 65536,
         dict_factor = 16,
         len_limit = 16,
         num_prev_pos = 1 << 16,
         pos_array_factor = 1 };

  int key4;			// key made from latest 4 bytes

public:
  explicit Fmatchfinder( const int ifd )
    :
    Matchfinder_base( before, dict_size, dict_factor, len_limit,
                      num_prev_pos, ifd, pos_array_factor ),
    key4( 0 )
    {}

  void reset() { Matchfinder_base::reset(); key4 = 0; }
  int longest_match_len( int * const distance );
  void longest_match_len();
  };


class FLZ_encoder : public LZ_encoder_base
  {
  Fmatchfinder & fmatchfinder;

  void move_pos( int n )
    {
    if( --n >= 0 ) fmatchfinder.move_pos();
    while( --n >= 0 )
      {
      fmatchfinder.longest_match_len();
      fmatchfinder.move_pos();
      }
    }

  void sequence_optimizer( int reps[num_rep_distances], State & state );

public:
  FLZ_encoder( Fmatchfinder & mf, const File_header & header, const int outfd )
    :
    LZ_encoder_base( header, mf.dictionary_size(), mf.match_len_limit(), outfd ),
    fmatchfinder( mf )
    {}

  bool encode_member( const long long member_size );
  };
