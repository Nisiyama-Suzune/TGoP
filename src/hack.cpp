/*	Hacks :
		unusual codes that may help in contests.
*/

//	long long formatting hack :

#ifdef WIN32
	#define LL "%I64d"
#else
	#define LL "%lld"
#endif

//	Optimizing hack :

#pragma GCC optimize ("O3")
#pragma GCC optimize ("whole-program")

//	Stack hack :

//	C++

#pragma comment(linker, "/STACK:36777216")

//	G++

int __size__ = 256 << 20;	//	256MB
char *__p__ = (char*)malloc(__size__) + __size__;
__asm__("movl %0, %%esp\n" :: "r"(__p__));

//	Ultra fast functions :

__inline void make_min (int &a, int &b) {
	asm (
		"cmpl %2, %0\n\t"
		"jle DONE\n\t"
		"movl %2, %0\n\t"
		"DONE:"
		: "=r" (a)
		: "0" (a), "r" (b)
	);
}

__inline void make_max (int &a, int &b) {
	asm (
		"cmpl %2, %0\n\t"
		"jge DONE\n\t"
		"movl %2, %0\n\t"
		"DONE:"
		: "=r" (a)
		: "0" (a), "r" (b)
	);
}

__inline int cmp (int a) {
	return (a >> 31) + (-a >> 31 & 1);
}

__inline int abs (int x) {
	int y = x >> 31;
	return (x + y) ^ y;
}

__inline int mul_mod (int a, int b) {
	int ret;
	asm (
		"mull %%ebx\n\t"
		"divl %%ecx\n\t"
		: "=d" (ret)
		: "a" (a), "b" (b), "c" (MO)
	);
	return ret;
}

