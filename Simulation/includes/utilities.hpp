#ifndef UTILITIES_H
#define UTILITIES_H

// Compute power of two greater than or equal to `n`
// https://www.techiedelight.com/round-next-highest-power-2/
// This only works for 32 bit integer
inline __attribute__ ((const)) uint findNextPowerOf2(uint n)
{
  // Might be able to remove this branch later
  if ( n==0 ) return 1;

  // decrement `n` (to handle the case when `n` itself is a power of 2)
  n--;

  // set all bits after the last set bit
  n |= n >> 1;
  n |= n >> 2;
  n |= n >> 4;
  n |= n >> 8;
  n |= n >> 16;

  // increment `n` and return
  return ++n;
}

#endif
