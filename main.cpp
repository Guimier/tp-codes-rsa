//
// TP6_RSA
//

#include <stdio.h>
#include <iostream>
#include <gmp.h>

#define BITSTRENGTH 14               /* size of modulus (n) in bits */
#define PRIMESIZE (BITSTRENGTH / 2)  /* size of the primes p and q  */

#define GN_POWM rew_powm
#define GN_NEXTPRIME rew_nextprime
#define GN_INVERT mpz_invert

#define SUSPICION 30

typedef std::pair<mpz_srcptr, mpz_srcptr> rsakey_t;
typedef unsigned long int ulong;

/*----- Rewritten functions -----*/

// Not managing the exp < 0 case since RSA does not need it.
void rew_powm( mpz_t rop, const mpz_t base, const mpz_t exp, const mpz_t mod )
{
	// rop = 1
	mpz_set_ui( rop, 1 );
	
	mpz_t rest;
	mpz_init( rest ); // portion of the exponent not yet applied
	// rest = exp
	mpz_set( rest, exp );
	
	mpz_t partiallyExponented;
	mpz_init( partiallyExponented );
	// partiallyExponented = base % mod
	mpz_mod( partiallyExponented, base, mod );
	
	mpz_t tmp;
	mpz_init( tmp ); // temporary var for intermediate values.
	
	// mpz_sgn is an easy way to test ( rest != 0 )
	while ( mpz_sgn( rest ) ) {
		if ( mpz_odd_p( rest ) ) {
			// tmp = rop * partiallyExponented
			mpz_mul( tmp, rop, partiallyExponented );
			// rop = tmp % mod
			mpz_mod( rop, tmp, mod );
		}
		// tmp = rest >> 2
		mpz_tdiv_q_2exp( rest, rest, 1 );
		// rest = tmp
		// mpz_set( rest, tmp );
		
		// tmp = partiallyExponented ^^ 2
		mpz_mul( tmp, partiallyExponented, partiallyExponented );
		// partiallyExponented = tmp % mod
		mpz_mod( partiallyExponented, tmp, mod );
	}
	
	mpz_clear( tmp );
	mpz_clear( partiallyExponented );
}

void loc_random( mpz_t rop, gmp_randstate_t state, ulong min, mpz_t max )
{
	mpz_t size;
	mpz_init( size );
	
	mpz_sub_ui( size, max, min );
	mpz_add_ui( size, size, 1 );
	
	mpz_urandomm( rop, state, size );
	mpz_add_ui( rop, rop, min );
	
	mpz_clear( size );
}

bool isCompositeRabinMiller( ulong s, const mpz_t t, const mpz_t a, const mpz_t nm1, const mpz_t n )
{
	mpz_t x;
	mpz_init( x );
	GN_POWM( x, a, t, n );
	bool res = true;
	
	if ( ! mpz_cmp_ui( x, 1 ) || ! mpz_cmp( x, nm1 ) ) {
		res = false;
		goto clearAndReturn;
	}
	for ( ulong r = 1; r < s; ++r ) {
		mpz_mul( x, x, x );
		mpz_mod( x, x, n );
		if ( ! mpz_cmp( x, nm1 ) ) {
			res = false;
			goto clearAndReturn;
		}
	}

clearAndReturn:
	mpz_clear( x );
	return res;
}

bool isProbablyPrimeRabinMiller( mpz_t n, ulong k )
{
	gmp_randstate_t randState;
	gmp_randinit_mt( randState );
	
	mpz_t nm1;
	mpz_init( nm1 );
	mpz_sub_ui( nm1, n, 1 ); // n minus 1
	
	mpz_t t;
	ulong s = 0;
	mpz_init( t );
	mpz_set( t, nm1 );
	// t is even
	while ( ! mpz_tstbit( t, 0 ) ) {
		++s;
		// t = t >> 2
		mpz_tdiv_q_2exp( t, t, 1 );
	}
	
	mpz_t a;
	mpz_init( a );
	
	bool res = true;
	ulong i = 0;
	while ( res && i < k ) {
		loc_random( a, randState, 2, nm1 );
		res = ! isCompositeRabinMiller( s, t, a, nm1, n );
		++i;
	}
	
	mpz_clear( a );
	mpz_clear( nm1 );
	mpz_clear( t );
	
	return res;
}

void rew_nextprime( mpz_t rop, const mpz_t init )
{
	mpz_set( rop, init );
	if ( mpz_even_p( rop ) ) {
		mpz_sub_ui( rop, rop, 1 );
	}
	if ( mpz_cmp_ui( rop, 1 ) < 0 ) {
		mpz_set_ui( rop, 1 );
	}
	
	do {
		mpz_add_ui( rop, rop, 2 );
	} while ( ! isProbablyPrimeRabinMiller( rop, SUSPICION ) );
}

void rew_invert (mpz_ptr res, mpz_srcptr b, mpz_srcptr n)
{
	mpz_t b0, n0, t0, t, q, r, tmp1, tmp2, tmp3, tmp4;
	
	mpz_init( b0 );
	mpz_init( n0 );
	mpz_init( t0 );
	mpz_init( t );
	mpz_init( q );
	mpz_init( r );
	mpz_init( tmp1 );
	mpz_init( tmp2 );
	mpz_init( tmp3 );
	mpz_init( tmp4 );
	
	mpz_set( b0,  b );
	mpz_set( n0,  n );
	
	mpz_set_ui( t,  1 );
	mpz_fdiv_q( q, n, b );
	mpz_mul( tmp1, q, b0 );
	mpz_sub( r, n0, tmp1 );
	
	while (mpz_cmp_si(r, 0) > 0)
	{	
		mpz_mul ( tmp2,  q,  t );
		mpz_sub( tmp2, t0, tmp2 );
			if( mpz_cmp_si( tmp2, 0) > 0 )
			{
				mpz_cdiv_r( tmp2, tmp2, n );
			}
			else
			{	
				mpz_neg( tmp2, tmp2 );
				mpz_fdiv_r( tmp3,  tmp2, n );
				mpz_sub( tmp2, n, tmp3 );
			}
		mpz_set( t0,  t);
		mpz_set( t, tmp2 );
		mpz_set( n0, b0 );
		mpz_set( b0, r );
		mpz_fdiv_q( q, n0, b0 );
		mpz_mul( tmp4, q, b0 );
		mpz_sub( r, n0, tmp4 );
	}

	if( mpz_cmp_si( t, 0) < 0 )
	{
		mpz_add(t, t, n);
	}
	mpz_set(res, t);
	
	mpz_clear( b0 );
	mpz_clear( n0 );
	mpz_clear( t0 );
	mpz_clear( t );
	mpz_clear( q );
	mpz_clear( r );
	mpz_clear( tmp1 );
	mpz_clear( tmp2 );
	mpz_clear( tmp3 );
	mpz_clear( tmp4 );
}

/*----- Helpers -----*/

std::ostream& operator << ( std::ostream& out, const mpz_t& nbr )
{
	return out << mpz_srcptr( nbr );
}

std::ostream& operator << ( std::ostream& out, mpz_srcptr nbr )
{
	size_t size = mpz_sizeinbase( nbr, 10 );
	char* str = new char[size+1];
	
	mpz_get_str( str, 10, nbr );
	out << str;
	
	delete[] str;
	
	return out;
}

std::ostream& operator << ( std::ostream& out, rsakey_t key)
{
	out << "( " << key.first << ", " << key.second << " )";
	return out;
}

/*----- RSA routines -----*/

void generate_e ( mpz_ptr e, gmp_randstate_t randState, mpz_srcptr phi )
{

	mpz_t gcd;
	mpz_init(gcd);
	
	do
	{
		mpz_urandomm( e, randState, phi);
		mpz_gcd( gcd, e, phi );
	
	} while ( mpz_cmp_ui( gcd, 1 ) );
	mpz_clear(gcd);
}


void encrypt (mpz_ptr res, mpz_srcptr message, rsakey_t key)
{
	GN_POWM(res, message, key.first, key.second);
}

#define decrypt encrypt

/*----- Main algorithme -----*/

>>>>>>> origin/master
/* Main subroutine */
int main( int, char** )
{

	mp_bitcnt_t bits = 64;
	mpz_t d, e, n;
	mpz_t tmp1, tmp2;
	mpz_t M, C, MM;
	
	/* Initialize the GMP integers */
	mpz_init( d );
	mpz_init( e );
	mpz_init( n );

	mpz_init( tmp1 );
	mpz_init( tmp2 );
	mpz_init( M );
	mpz_init( C );
	mpz_init( MM );
	

	gmp_randstate_t randState;
	gmp_randinit_mt( randState );

	/* This function creates the keys. The basic algorithm is...
	 *
	 *  1. Generate two large distinct primes p and q randomly
	 *  2. Calculate n = pq and x = (p-1)(q-1)
	 *  3. Select a random integer e (1<e<x) such that gcd(e, x) = 1
	 *  4. Calculate the unique d such that ed = 1(mod x)
	 *  5. Public key pair : (e, n), Private key pair : (d, n)
	 *
	 */

	/*
	 *  Step 1 : Get two large primes.
	 */
	mpz_t p, q;
	mpz_init( p );
	mpz_init( q );
	
	mpz_urandomb( tmp1, randState, bits );
	GN_NEXTPRIME( p, tmp1 );
	mpz_urandomb( tmp1, randState, bits );
	GN_NEXTPRIME( q, tmp1 );

	std::cout << "Random Prime 'p' = " << p <<  std::endl;
	std::cout << "Random Prime 'q' = " << q <<  std::endl;

	/*
	 *  Step 2 : Calculate n (=pq) ie the 1024 bit modulus
	 *  and x (=(p-1)(q-1)).
	 */
	mpz_t phi;
	mpz_init( phi );

	/* Calculate n... */
	mpz_mul( n, p, q );
	std::cout << "\t n = " << n << std::endl;


	/* Calculate x... */
	mpz_t p_minus_1, q_minus_1;
	mpz_init( p_minus_1 );
	mpz_init( q_minus_1 );

	mpz_sub_ui( p_minus_1, p, (unsigned long int)1 );
	mpz_sub_ui( q_minus_1, q, (unsigned long int)1 );

	mpz_mul( phi, p_minus_1, q_minus_1 );
	std::cout << "\t phi(n) = " << phi << std::endl;

	/*
	 *  Step 3 : Get small odd integer e such that gcd(e, x) = 1.
	 */
	//mpz_init_set_str( e, "79", 0 );
	generate_e( e, randState, phi );
	std::cout << "\t e = " << e << std::endl;

	/*
	 *  Step 4 : Calculate unique d such that ed = 1(mod x)
	 */
	//mpz_init_set_str( d, "1019", 0 );
	GN_INVERT( d, e, phi );
	std::cout << "\t d = " << d << std::endl << std::endl;
	
	rsakey_t pub = std::make_pair( e, n );
	rsakey_t priv = std::make_pair( d, n );
	
	/*
	 *  Step 5 : Print the public and private key pairs...
	 */
	std::cout << "Public Keys  (e, n): " << pub << std::endl;
	std::cout << "Private Keys (d, n): " << priv << std::endl << std::endl;
	
	/*
	 *  Encrypt
	 */
	mpz_urandomm( M, randState, n );
	std::cout << "Message: " << M << std::endl;
	encrypt(C, M, pub);
	std::cout << "Encrypted: "<< C << std::endl;
	decrypt( MM, C, priv );
	std::cout << "Decrypted: "<< MM << std::endl;
	
	/* Clean up the GMP integers */
	mpz_clear( p_minus_1 );
	mpz_clear( q_minus_1 );
	mpz_clear( phi );
	mpz_clear( p );
	mpz_clear( q );

	mpz_clear( d );
	mpz_clear( e );
	mpz_clear( n );

	mpz_clear( tmp1 );
	mpz_clear( tmp2 );
	mpz_clear( M );
	mpz_clear( C );
	mpz_clear( MM );
	
	return 0;
}

