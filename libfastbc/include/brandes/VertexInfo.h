#ifndef FASTBC_BRANDES_VERTEXINFO_H
#define FASTBC_BRANDES_VERTEXINFO_H

#include <algorithm>
#include <cmath>
#include <vector>

#ifndef FASTBC_BRANDES_VERTEXINFO_PENALTY
#define FASTBC_BRANDES_VERTEXINFO_PENALTY 1000
#endif

namespace fastbc {
	namespace brandes {

		template<typename V, typename W>
		class VertexInfo
		{
		public:
			/**
			 *	@brief Initialize a Vertex info object with borderCount borders info
			 * 
			 *	@param borderCount Number of border info to store
			 */
			VertexInfo(int borderCount);

			template<typename N, typename E>
			VertexInfo(const VertexInfo<N, E>& copy);

			template<typename N, typename E>
			VertexInfo<V, W>& operator=(const VertexInfo<N, E>& other);

			void setBorderSPLength(int storeIndex, W length);

			W getBorderSPLength(int storeIndex) const;

			void setBorderSPCount(int storeIndex, V count);

			V getBorderSPCount(int storeIndex) const;

			W getMinBorderSPLength() const;

			void normalize();

			void reset();

			int borders() const;

			template<typename N, typename E>
			W squaredDistance(const VertexInfo<N, E>& other) const;

			template<typename N, typename E>
			W contributionDistance(const VertexInfo<N, E>& other) const;

			template<typename N, typename E>
			VertexInfo<V, W>& operator+=(const VertexInfo<N, E>& other);

			template<typename N, typename E>
			VertexInfo<V, W>& operator-=(const VertexInfo<N, E>& other);

			template<typename N, typename E>
			VertexInfo<V, W>& operator*=(const VertexInfo<N, E>& other);

			template<typename N, typename E>
			VertexInfo<V, W>& operator/=(const VertexInfo<N, E>& other);

			template<typename N, typename E>
			VertexInfo<V, W> operator+(const VertexInfo<N, E>& other);

			template<typename N, typename E>
			VertexInfo<V, W> operator-(const VertexInfo<N, E>& other);

			template<typename N, typename E>
			VertexInfo<V, W> operator*(const VertexInfo<N, E>& other);

			template<typename N, typename E>
			VertexInfo<V, W> operator/(const VertexInfo<N, E>& other);

			template<typename T> VertexInfo<V, W>& operator+=(T num);

			template<typename T> VertexInfo<V, W>& operator-=(T num);

			template<typename T> VertexInfo<V, W>& operator*=(T num);

			template<typename T> VertexInfo<V, W>& operator/=(T num);

			template<typename T> VertexInfo<V, W> operator+(T num);

			template<typename T> VertexInfo<V, W> operator-(T num);

			template<typename T> VertexInfo<V, W> operator*(T num);

			template<typename T> VertexInfo<V, W> operator/(T num);

			template<typename N, typename E>
			bool operator==(const VertexInfo<N, E>& other) const;

			template<typename N, typename E>
			bool operator!=(const VertexInfo<N, E>& other) const;

			template<typename N, typename E>
			bool operator<(const VertexInfo<N, E>& other) const;

			template<typename N, typename E>
			bool operator>(const VertexInfo<N, E>& other) const;

			template<typename N, typename E>
			bool operator<=(const VertexInfo<N, E>& other) const;

			template<typename N, typename E>
			bool operator>=(const VertexInfo<N, E>& other) const;

			template<typename N, typename E>
			W compare(const VertexInfo<N, E>& other) const;

		private:
			int _borderCount;
			std::vector<W> _borderSPLength;
			std::vector<V> _borderSPCount;

			template<typename N, typename E>
			friend class VertexInfo;
		};

	}
}

template<typename V, typename W>
fastbc::brandes::VertexInfo<V, W>::VertexInfo(int borderCount)
	: _borderCount(borderCount),
	_borderSPLength(_borderCount),
	_borderSPCount(_borderCount)
{
}

template<typename V, typename W>
template<typename N, typename E>
fastbc::brandes::VertexInfo<V, W>::VertexInfo(const VertexInfo<N, E>& copy)
	: _borderCount(copy._borderCount),
	_borderSPLength(_borderCount),
	_borderSPCount(_borderCount)
{
	#pragma omp simd
	for (int i = 0; i < _borderCount; ++i)
	{
		_borderSPLength[i] = (W)copy._borderSPLength[i];
		_borderSPCount[i] = (V)copy._borderSPCount[i];
	}
}

template<typename V, typename W>
template<typename N, typename E>
fastbc::brandes::VertexInfo<V, W>& fastbc::brandes::VertexInfo<V, W>::operator=(const VertexInfo<N, E>& other)
{
	if (_borderCount != other._borderCount)
	{
		_borderCount = other._borderCount;
		_borderSPLength.resize(_borderCount);
		_borderSPCount.resize(_borderCount);
	}

	#pragma omp simd
	for (int i = 0; i < _borderCount; ++i)
	{
		_borderSPLength[i] = (W)other._borderSPLength[i];
		_borderSPCount[i] = (V)other._borderSPCount[i];
	}

	return *this;
}

template<typename V, typename W>
void fastbc::brandes::VertexInfo<V, W>::setBorderSPLength(int storeIndex, W length)
{
	if (storeIndex < _borderCount)
	{
		_borderSPLength[storeIndex] = length;
	}
	else
	{
		throw std::out_of_range("Given store index is out of range.");
	}
}

template<typename V, typename W>
W fastbc::brandes::VertexInfo<V, W>::getBorderSPLength(int storeIndex) const
{
	if (storeIndex < _borderCount)
	{
		return _borderSPLength[storeIndex];
	}
	else
	{
		throw std::out_of_range("Given store index is out of range.");
	}
}

template<typename V, typename W>
void fastbc::brandes::VertexInfo<V, W>::setBorderSPCount(int storeIndex, V count)
{
	if (storeIndex < _borderCount)
	{
		_borderSPCount[storeIndex] = count;
	}
	else
	{
		throw std::out_of_range("Given store index is out of range.");
	}
}

template<typename V, typename W>
V fastbc::brandes::VertexInfo<V, W>::getBorderSPCount(int storeIndex) const
{
	if (storeIndex < _borderCount)
	{
		return _borderSPCount[storeIndex];
	}
	else
	{
	throw std::out_of_range("Given store index is out of range.");
	}
}

template<typename V, typename W>
W fastbc::brandes::VertexInfo<V, W>::getMinBorderSPLength() const
{
	// It could be possible to have a sub-grph not connected to external vertices
	if (!_borderCount) { return 0; }

	return *std::min_element(_borderSPLength.begin(), _borderSPLength.end());
}

template<typename V, typename W>
void fastbc::brandes::VertexInfo<V, W>::normalize()
{
	W min = getMinBorderSPLength();

	#pragma omp simd
	for (int i = 0; i < _borderCount; ++i)
	{
		_borderSPLength[i] -= min;
	}
}

template<typename V, typename W>
void fastbc::brandes::VertexInfo<V, W>::reset()
{
	_borderSPLength.assign(_borderCount, (W)0);
	_borderSPCount.assign(_borderCount, (V)0);
}

template<typename V, typename W>
int fastbc::brandes::VertexInfo<V, W>::borders() const
{
	return _borderCount;
}

template<typename V, typename W>
template<typename N, typename E>
W fastbc::brandes::VertexInfo<V, W>::squaredDistance(const VertexInfo<N, E>& other) const
{
	W sqDistance = 0;

	#pragma omp simd reduction(+:sqDistance)
	for (int i = 0; i < _borderCount; ++i)
	{
		sqDistance += std::pow(_borderSPLength[i] - (W)other._borderSPLength[i], 2);
		sqDistance += std::pow(_borderSPCount[i] - (V)other._borderSPCount[i], 2);
	}

	return sqDistance;
}

template<typename V, typename W>
template<typename N, typename E>
W fastbc::brandes::VertexInfo<V, W>::contributionDistance(const VertexInfo<N, E>& other) const
{
	W cDistance = 0;

	#pragma omp simd reduction(+:cDistance)
	for (int i = 0; i < _borderCount; ++i)
	{
		if (_borderSPCount[i] != 0 || other._borderSPCount[i] != 0)
		{
			if (_borderSPCount[i] > 0 && other._borderSPCount[i] > 0)
			{
				cDistance += std::pow(_borderSPLength[i] - (W)other._borderSPLength[i], 2);
				cDistance += std::pow(_borderSPCount[i] - (V)other._borderSPCount[i], 2);
			}
			else
			{
				cDistance += (W)FASTBC_BRANDES_VERTEXINFO_PENALTY;
			}
		}
	}

	return cDistance;
}

template<typename V, typename W>
template<typename N, typename E>
fastbc::brandes::VertexInfo<V, W>&
fastbc::brandes::VertexInfo<V, W>::operator+=(const VertexInfo<N, E>& other)
{
	#pragma omp simd
	for (int i = 0; i < _borderCount; ++i)
	{
		_borderSPLength[i] += other._borderSPLength[i];
		_borderSPCount[i] += other._borderSPCount[i];
	}

	return *this;
}

template<typename V, typename W>
template<typename N, typename E>
fastbc::brandes::VertexInfo<V, W>& 
fastbc::brandes::VertexInfo<V, W>::operator-=(const VertexInfo<N, E>& other)
{
	#pragma omp simd
	for (int i = 0; i < _borderCount; ++i)
	{
		_borderSPLength[i] -= other._borderSPLength[i];
		_borderSPCount[i] -= other._borderSPCount[i];
	}

	return *this;
}

template<typename V, typename W>
template<typename N, typename E>
fastbc::brandes::VertexInfo<V, W>&
fastbc::brandes::VertexInfo<V, W>::operator*=(const VertexInfo<N, E>& other)
{
	#pragma omp simd
	for (int i = 0; i < _borderCount; ++i)
	{
		_borderSPLength[i] *= other._borderSPLength[i];
		_borderSPCount[i] *= other._borderSPCount[i];
	}

	return *this;
}

template<typename V, typename W>
template<typename N, typename E>
fastbc::brandes::VertexInfo<V, W>&
fastbc::brandes::VertexInfo<V, W>::operator/=(const VertexInfo<N, E>& other)
{
	#pragma omp simd
	for (int i = 0; i < _borderCount; ++i)
	{
		_borderSPLength[i] /= other._borderSPLength[i];
		_borderSPCount[i] /= other._borderSPCount[i];
	}

	return *this;
}

template<typename V, typename W>
template<typename N, typename E>
fastbc::brandes::VertexInfo<V, W>
fastbc::brandes::VertexInfo<V, W>::operator+(const VertexInfo<N, E>& other)
{
	VertexInfo<V, W> sum(*this);

	return sum += other;
}

template<typename V, typename W>
template<typename N, typename E>
fastbc::brandes::VertexInfo<V, W> 
fastbc::brandes::VertexInfo<V, W>::operator-(const VertexInfo<N, E>& other)
{
	VertexInfo<V, W> sub(*this);

	return sub -= other;
}

template<typename V, typename W>
template<typename N, typename E>
fastbc::brandes::VertexInfo<V, W>
fastbc::brandes::VertexInfo<V, W>::operator*(const VertexInfo<N, E>& other)
{
	VertexInfo<V, W> mul(*this);

	return mul *= other;
}

template<typename V, typename W>
template<typename N, typename E>
fastbc::brandes::VertexInfo<V, W>
fastbc::brandes::VertexInfo<V, W>::operator/(const VertexInfo<N, E>& other)
{
	VertexInfo<V, W> div(*this);

	return div /= other;
}

template<typename V, typename W>
template<typename T>
fastbc::brandes::VertexInfo<V, W>&
fastbc::brandes::VertexInfo<V, W>::operator+=(T num)
{
	#pragma omp simd
	for (int i = 0; i < _borderCount; ++i)
	{
		_borderSPLength[i] += num;
		_borderSPCount[i] += num;
	}

	return *this;
}

template<typename V, typename W>
template<typename T>
fastbc::brandes::VertexInfo<V, W>&
fastbc::brandes::VertexInfo<V, W>::operator-=(T num)
{
	#pragma omp simd
	for (int i = 0; i < _borderCount; ++i)
	{
		_borderSPLength[i] -= num;
		_borderSPCount[i] -= num;
	}

	return *this;
}

template<typename V, typename W>
template<typename T>
fastbc::brandes::VertexInfo<V, W>&
fastbc::brandes::VertexInfo<V, W>::operator*=(T num)
{
	#pragma omp simd
	for (V i = 0; i < _borderCount; ++i)
	{
		_borderSPLength[i] *= num;
		_borderSPCount[i] *= num;
	}

	return *this;
}

template<typename V, typename W>
template<typename T>
fastbc::brandes::VertexInfo<V, W>&
fastbc::brandes::VertexInfo<V, W>::operator/=(T num)
{
	#pragma omp simd
	for (int i = 0; i < _borderCount; ++i)
	{
		_borderSPLength[i] /= num;
		_borderSPCount[i] /= num;
	}

	return *this;
}

template<typename V, typename W>
template<typename T>
fastbc::brandes::VertexInfo<V, W>
fastbc::brandes::VertexInfo<V, W>::operator+(T num)
{
	VertexInfo<V, W> sum(*this);

	return sum += num;
}

template<typename V, typename W>
template<typename T>
fastbc::brandes::VertexInfo<V, W>
fastbc::brandes::VertexInfo<V, W>::operator-(T num)
{
	VertexInfo<V, W> sub(*this);

	return sub -= num;
}

template<typename V, typename W>
template<typename T>
fastbc::brandes::VertexInfo<V, W>
fastbc::brandes::VertexInfo<V, W>::operator*(T num)
{
	VertexInfo<V, W> mul(*this);

	return mul *= num;
}

template<typename V, typename W>
template<typename T>
fastbc::brandes::VertexInfo<V, W>
fastbc::brandes::VertexInfo<V, W>::operator/(T num)
{
	VertexInfo<V, W> div(*this);

	return div /= num;
}

template<typename V, typename W>
template<typename N, typename E>
bool fastbc::brandes::VertexInfo<V, W>::operator==(const VertexInfo<N, E>& other) const
{
	return compare(other) == 0;
}

template<typename V, typename W>
template<typename N, typename E>
bool fastbc::brandes::VertexInfo<V, W>::operator!=(const VertexInfo<N, E>& other) const
{
	return compare(other) != 0;
}

template<typename V, typename W>
template<typename N, typename E>
bool fastbc::brandes::VertexInfo<V, W>::operator<(const VertexInfo<N, E>& other) const
{
	return compare(other) < 0;
}

template<typename V, typename W>
template<typename N, typename E>
bool fastbc::brandes::VertexInfo<V, W>::operator>(const VertexInfo<N, E>& other) const
{
	return compare(other) > 0;
}

template<typename V, typename W>
template<typename N, typename E>
bool fastbc::brandes::VertexInfo<V, W>::operator<=(const VertexInfo<N, E>& other) const
{
	return compare(other) <= 0;
}

template<typename V, typename W>
template<typename N, typename E>
bool fastbc::brandes::VertexInfo<V, W>::operator>=(const VertexInfo<N, E>& other) const
{
	return compare(other) >= 0;
}

template<typename V, typename W>
template<typename N, typename E>
W fastbc::brandes::VertexInfo<V, W>::compare(const VertexInfo<N, E>& other) const
{
	for (int i = 0; i < _borderCount; i++)
	{
		if (W cmp = _borderSPCount[i] - other._borderSPCount[i]; cmp != 0)
		{
			return cmp;
		}

		if (W cmp = _borderSPLength[i] - other._borderSPLength[i]; cmp != 0)
		{
			return cmp;
		}
	}

	return 0;
}

#endif // !FASTBC_BRANDES_IVERTEXINFO_H
