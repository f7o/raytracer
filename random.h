#ifndef RANDOM_H
#define RANDOM_H

// George Marsaglia based pseudo-random number generator
// faster than std::rand()
class RandomGenerator {

private:
	// state variables
	unsigned int x, y, z, w;

public:
	void initialize(unsigned int seed) {
		x = (123456789u ^ seed) * 88675123u;
		y = (362436069u ^ seed) * 123456789u;
		z = (521288629u ^ seed) * 362436069u;
		w = (88675123u ^ seed) * 521288629u;

		randInt();
	}

	void initialize(unsigned int a,
	                unsigned int b,
	                unsigned int c,
	                unsigned int d) {
		x = a;
		y = b;
		z = c;
		w = d;

		randInt();
	}

	unsigned int randInt() {
		unsigned int t;

		t = x ^ (x << 11u);
		x = y;
		y = z;
		z = w;
		return w = w ^ (w >> 19u) ^ (t ^ (t >> 8u));
	}

	float randFloat() {
		return 2.32830643653869629E-10f * randInt();
	}
};

#endif // RANDOM_H
