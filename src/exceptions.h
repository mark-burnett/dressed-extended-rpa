/* Exceptions used throughout library.
 *
 * Mark Burnett, November 2008
 */

#ifndef _my_exceptions_h_
#define _my_exceptions_h_

#include <stdexcept>

class particle_number_mismatch : public std::exception {};

class illegal_quantum_number   : public std::exception { };

class illegal_angular_momentum : public illegal_quantum_number { };
class illegal_parity           : public illegal_quantum_number { };
class illegal_isospin          : public illegal_quantum_number { };

class map_key_missing : public std::exception { };

class file_error : public std::exception { };

class bad_comparison : public std::exception { };

class root_finding_error : public std::exception { };

class invalid_matrix_position : public std::exception { };

#endif // _my_exceptions_h_
