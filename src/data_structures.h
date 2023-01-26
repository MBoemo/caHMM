//----------------------------------------------------------
// Copyright 2023 University of Cambridge
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under the MIT license.
// You should have received a copy of the license with
// this software.  If not, please Email the author.
//----------------------------------------------------------

#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

#include <stdint.h>
#include <string.h> /*memset */
#include "error_handling.h"


/*this is gross, but fix later */
/*this is gross, but fix later */
class sparseMatrix {

	private:
		size_t _rows, _cols;
		std::map< std::pair< size_t, size_t >, double > _data;

	public:
		sparseMatrix( size_t rows, size_t cols ){

			_rows = rows;
			_cols = cols;
		}
		~sparseMatrix(){
			_data.clear();
		}

		inline void assign( size_t row, size_t col, double value ){

			if ( row >= _rows or col >= _cols ) throw BadIndex();
			_data[ std::make_pair(row,col) ] = value;
		}

		inline double access( size_t row, size_t col ){

			if ( row >= _rows or col >= _cols ) throw BadIndex();
			
			std::pair< size_t, size_t > key = std::make_pair( row, col );

			auto i = _data.find( key );
			if ( i == _data.end() ) return NAN;
			else return _data[ key ];
		}
};




class defaultSparse {

	private:
		size_t _rows, _cols;
		std::map< std::pair< size_t, size_t >, double > _data;

	public:
		defaultSparse( size_t rows, size_t cols ){

			_rows = rows;
			_cols = cols;
		}
		~defaultSparse(){
			_data.clear();
		}

		inline void assign( size_t row, size_t col, double value ){

			if ( row >= _rows or col >= _cols ) throw BadIndex();
			_data[ std::make_pair(row,col) ] = value;
		}

		inline double access( size_t row, size_t col ){

			if ( row >= _rows or col >= _cols ) throw BadIndex();
			
			std::pair< size_t, size_t > key = std::make_pair( row, col );

			auto i = _data.find( key );
			if ( i == _data.end() ) return NAN;
			else return _data[ key ];
		}
};

/*struct for a matrix-like object that uses a contiguous block of memory */
template< typename type >
struct Matrix {
	type *cells;
	uint32_t numRows;
	uint32_t numCols;
};


/*for now, restrict attention to matrices of doubles and ints */
typedef Matrix< double > DoubleMatrix;
typedef Matrix< int > IntMatrix;


/*allocate a contiguous block of memory for the matrix */
template< typename type >
void allocate( Matrix< type >& matrix, uint32_t rows, uint32_t cols )
{
	matrix.numRows = rows;
	matrix.numRows = cols;
    
	uint32_t N = matrix.numRows * matrix.numCols;
	matrix.cells = ( type* ) malloc( N * sizeof( type ) );
	memset( matrix.cells, 0, N * sizeof( type ) );
}


/*return the appropriate cell, given row and column indicies */
template< typename type > 
inline uint32_t cell( const Matrix< type >& matrix, uint32_t row, uint32_t col )
{
	return row * matrix.numCols + col;
}


/*set a cell in the matrix to a value, using its row and column indices */
template< typename type, typename U >
inline void set( Matrix< type >& matrix, uint32_t row, uint32_t col, U value )
{
	uint32_t c = cell( matrix, row, col ); /*grab the cell that we're interested in */
	matrix.cells[ c ] = value; /*set that cell to the value that we want */
}


/*retrive a value from the matrix at a certain cell position, given by row and column indicies */
template< typename type >
inline type get( const Matrix< type >& matrix, uint32_t row, uint32_t col )
{
	uint32_t c = cell( matrix, row, col );
	return matrix.cells[ c ];
}

#endif
