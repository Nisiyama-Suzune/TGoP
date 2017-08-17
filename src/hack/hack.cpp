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

#define __ __attribute__ ((optimize ("-O3"))) 
#define _ __ __inline __attribute__ ((__gnu_inline__, __always_inline__, __artificial__)) __attribute__ ((aligned))
#define NDEBUG

//	Stack hack :

//	C++

#pragma comment(linker, "/STACK:36777216")

//	G++

int __size__ = 256 << 20;	//	256MB
char *__p__ = (char*)malloc(__size__) + __size__;
__asm__ ("movl %0, %%esp\n" :: "r"(__p__));

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

__inline int next_uint () {
	const int SIZE = 110000; static char buf[SIZE]; static int p = SIZE;
	register int ans = 0, f = 1;
	while ((p < SIZE || fread (buf, 1, SIZE, stdin) && (p = 0, 1))
		&& (isdigit (buf[p]) && (ans = ans * 10 + buf[p] - '0', f = 0, 1) || f)) ++p;
	return ans;
}

__inline int next_int () {
	const int SIZE = 110000; static char buf[SIZE]; static int p = SIZE;
	register int ans = 0, f = 1, sgn = 1;
	while ((p < SIZE || fread (buf, 1, SIZE, stdin) && (p = 0, 1)) && 
		(isdigit (buf[p]) && (ans = ans * 10 + buf[p] - '0', f = 0, 1) || 
		f && (buf[p] == '-' && (sgn = 0), 1))) ++p;
	return sgn ? ans : -ans;
}
