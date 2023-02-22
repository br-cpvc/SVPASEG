#if! defined(__RFRANDOM_H__)
#define __RFRANDOM_H__

#define RANDOM_N              (624>>1)              // Length of state vector
#define RANDOM_M              (397>>1)              // A period parameter
#define RANDOM_K              (0x9908B0DFU)         // A magic constant

class RFRandom
{
private:
	unsigned int	m_auiState[RANDOM_N + 1];	// State vector + 1 extra to not violate ANSI C
	unsigned int*	m_uiNext;					// Next random value is computed from here
	int		m_iLeft;					// Can *next++ this many times before reloading

	bool	m_bRegenerate;
	int		m_uiLast;			// TEP: Modification.. used for next seed!

	unsigned int Reload();

public:
	RFRandom();
	virtual ~RFRandom();

	void Seed(unsigned int uiSeed);

	// Get random 32bit integer
	unsigned int Get();

	// Get random 32bit float [-1.0f;1.0f[
	float GetFloatNorm();

	// Get random 32bit float [0.0f;1.0f[
	float GetFloatNormPositive();

};

#endif // __KLRANDOM_H__
