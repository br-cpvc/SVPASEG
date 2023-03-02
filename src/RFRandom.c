// This is the ``Mersenne Twister'' random number generator MT19937, which
// generates pseudorandom integers uniformly distributed in 0..(2^32 - 1)
// starting from any odd seed in 0..(2^32 - 1).  This version is a recode
// by Rafael Baptista, based on a recode by Shawn Cokus 
// (Cokus@math.washington.edu) on March 8, 1998 of a version by
// Takuji Nishimura (who had suggestions from Topher Cooper and Marc Rieffel in
// July-August 1997).
//
// According to the URL 
// (and paraphrasing a bit in places), the Mersenne Twister is ``designed
// with consideration of the flaws of various existing generators,'' has
// a period of 2^19937 - 1, gives a sequence that is 623-dimensionally
// equidistributed, and ``has passed many stringent tests, including the
// die-hard test of G. Marsaglia and the load test of P. Hellekalek and
// S. Wegenkittl.''  It is efficient in memory usage (typically using 2506
// to 5012 bytes of static data, depending on data type sizes, and the code
// is quite short as well).  It generates random numbers in batches of 624
// at a time, so the caching and pipelining of modern systems is exploited.
// It is also divide- and mod-free.
///
//
// The code as Shawn received it included the following notice:
//
//   Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.  When
//   you use this, send an e-mail to  with
//   an appropriate reference to your work.
//
// It would be nice to CC:  when you write.

#include "RFRandom.h"

#define hiBit(u)       ((u) & 0x80000000U)   // Mask all but highest   bit of u
#define loBit(u)       ((u) & 0x00000001U)   // Mask all but lowest    bit of u
#define loBits(u)      ((u) & 0x7FFFFFFFU)   // Mask     the highest   bit of u
#define mixBits(u, v)  (hiBit(u)|loBits(v))  // Move hi bit of u to hi bit of v

//-------------------------------------------------------------------------

RFRandom::RFRandom() :
	m_bRegenerate(true)
{
	Seed(0);
}

RFRandom::~RFRandom()
{
}

//-------------------------------------------------------------------------

void RFRandom::Seed(unsigned int uiSeed)
{
	unsigned int uiX	= (uiSeed | 0x1 ) & 0xffffffff;
	unsigned int* puiS	= m_auiState;
	unsigned int uiJ 	= RANDOM_N;
	m_iLeft		= 0;
	*puiS++		= uiX;

	while (--uiJ)
	{
		uiX 	*= 69069;
		*puiS = (uiX & 0xffffffff);
		puiS++;
	}

	m_bRegenerate = true;
}

unsigned int RFRandom::Reload()
{
	unsigned int* puiP0 = m_auiState;
	unsigned int* puiP2 = m_auiState + 2;
	unsigned int* puiPM = m_auiState + RANDOM_M;
	unsigned int uiS0;
	unsigned int uiS1;

	int iJ;

	// if this is the first time through seed the algorithm.
	if ( m_iLeft < -1 ) Seed(m_uiLast);

	m_iLeft = RANDOM_N-1;
	m_uiNext = &(m_auiState[1]);

	for( uiS0 = m_auiState[0], uiS1 = m_auiState[1], iJ = (RANDOM_N-RANDOM_M+1); --iJ; uiS0 = uiS1, uiS1 = *puiP2++ )
		*puiP0++ = ( *puiPM++ ^ (mixBits(uiS0, uiS1) >> 1) ^ (loBit(uiS1) ? RANDOM_K : 0));

	for( puiPM = m_auiState, iJ = RANDOM_M; --iJ; uiS0 = uiS1, uiS1 = *puiP2++ )
		*puiP0++ = ( *puiPM++ ^ ( mixBits(uiS0, uiS1) >> 1 ) ^ ( loBit(uiS1) ? RANDOM_K : 0 ));

	uiS1 = m_auiState[0];

	*puiP0 = ( *puiPM ^ ( mixBits( uiS0, uiS1 ) >> 1 ) ^ ( loBit( uiS1 ) ? RANDOM_K : 0 ));

	uiS1 ^= ( uiS1 >> 11 );
	uiS1 ^= ( uiS1 <<  7 ) & 0x9D2C5680U;
	uiS1 ^= ( uiS1 << 15 ) & 0xEFC60000U;

	return m_uiLast = ( uiS1 ^ ( uiS1 >> 18 ));
}

unsigned int RFRandom::Get()
{
	unsigned int uiY;

	m_iLeft--;
	if ( m_iLeft < 0 ) 
	{
		if (m_bRegenerate)
		{
			return Reload();
		} 
		else 
		{
			m_iLeft = RANDOM_N-1;
			m_uiNext = &(m_auiState[1]);
		}
	}

	uiY  = *m_uiNext;
	m_uiNext++;
	uiY ^= ( uiY >> 11 );
	uiY ^= ( uiY <<  7 ) & 0x9D2C5680;
	uiY ^= ( uiY << 15 ) & 0xEFC60000;
	uiY ^= ( uiY >> 18 );

	return m_uiLast = uiY;
}

float RFRandom::GetFloatNorm()
{
	float fVal = (float)(((int)Get() & 0x7fffffff) - 0x40000000);
	return fVal / ((float)(0x40000000));
}

float RFRandom::GetFloatNormPositive()
{
	float fVal = (float)(Get() & 0x7fffffff);
	return fVal / ((float)(0x80000000));
}

//-------------------------------------------------------------------------
