#ifndef TIMER_CLASS
#define TIMER_CLASS

/*
  Simple timer class from learncpp.com
*/

#include <chrono> // for std::chrono functions
#include <iostream>

class Timer
{
private:
	// Type aliases to make accessing nested type easier
	using clock_t = std::chrono::high_resolution_clock;
	using second_t = std::chrono::duration<double, std::ratio<1> >;

	std::chrono::time_point<clock_t> m_beg;
	double m_accumulate { 0.0 };

public:
	Timer() : m_beg(clock_t::now())
	{
	}

	void reset()
	{
		m_beg = clock_t::now();
	}

	double elapsed() const
	{
		return static_cast<double>(
			std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count()
		);
	}

	void accumulate()
	{
		m_accumulate += this->elapsed();
	}

	friend std::ostream& operator<< (std::ostream &out, const Timer &t)
	{
		out << t.m_accumulate;
		return out;
	}

	double getTotal() const
	{
		return m_accumulate;
	}
};

#endif
