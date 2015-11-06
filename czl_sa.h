#ifndef CZL_SA_H
#define CZL_SA_H

#include "czl_common.h"

/**
 * Ref. 
 * Nong, Ge; Zhang, Sen; Chan, Wai Hong (2009). Linear Suffix Array Construction 
 * by Almost Pure Induced-Sorting. 2009 Data Compression Conference. p. 193. 
 * doi:10.1109/DCC.2009.42. ISBN 978-0-7695-3592-0.
 * 
 * @param  n   size of text array
 * @param  t   content of text array
 * @param  sa  suffix array
 */
template<typename T, typename S=int>
int sa_is_Nong(S n, T* t, T* sa)
{
    bool* sl[n];  // false: S-type, true: T-type
	S LMS_n=0;
	T* LMS[n];
    /// Scan 't' once to classify all the characters as L- or S-type into 'sl';
	sa_is_Nong_sl(n, t, bool *sl);
	/// Scan 't' once to find all the LMS-substrings in S into LMS ;
	for (S i=1; i<n; i++) {
		if (sl[i-1] && !sl[i]) {
			LMS[LMS_n++] = t+i;
		}
	}
	/// Induced sort all the LMS-substrings using P1 and B;
	/// Name each LMS-substring in S by its bucket index to get a new shortened string t1;
	/// if Each character in t1 is unique
	/// then
	/// Directly compute SA1 from S1 ;
	/// else
	/// SA - IS(S 1 , SA 1 ); ✄ where recursive call happens
	/// Induce SA from SA 1 ;
	return 0;
}

template<typename T, typename S=int>
int sa_is_Nong_SL(S n, T* t, T* sa, bool *sl)
{
	if (n==0) return 0;
	S i=n;
	sl[i] = false;
	sl[--i] = true;
	while (i>0) {
		if (t[i-1] < t[i]) {
			sl[i-1] = false;
		} else if (t[i-1] == t[i]) {
			sl[i-1] = sl[i];
		} else {
			sl[i-1] = true;
		}
		i--;
	}
	return 0;
}
#endif
