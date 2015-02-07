#ifndef VECTOR_H
#define VECTOR_H

#include <math.h>
#include <ostream>
#include "matrix.h"

using namespace std;

template<typename Type, size_t ROW_NUM, size_t COL_NUM> 
class Vector : public Matrix<Type, ROW_NUM, COL_NUM>
{
public:
	template <typename... T>
	Vector(T... ts) : Matrix({ ts... }) { }

	Type getLengthSq() const { 
		Type sum;
		for (int i = 0; i < getElemNum(); i++){
			sum += elem(i)*elem(i);
		}
		
		return sum;
	}

	operator Matrix<Type, ROW_NUM, COL_NUM>(){
		return Matrix<Type, ROW_NUM, COL_NUM>(data._Elems);
	}

	Type getLength(void) const{ 
		return (Type)sqrt(getLengthSq()); 
	}
	
	void normalize()
	{
		Type length = getLength();
		if (length == (Type)0) return;
		(*this) /= length;
	}

	Vector getNormalized() const
	{
		Vector<Type, ROW_NUM, COL_NUM> temp(*this);
		temp.normalize();
		return temp;
	}
};

template<typename Type2, size_t ROW_NUM, size_t COL_NUM>
Type2 operator*(const Vector<Type2, ROW_NUM, COL_NUM> &lhs, const Vector<Type2, ROW_NUM, COL_NUM> &rhs) {
	Type2 dotP = 0;
	for (int i = 0; i < lhs.getElemNum(); i++){
		dotP += lhs(i)*rhs(i);
	}

	return dotP;
}

template <typename Type2, typename Type, size_t ROW_NUM, size_t COL_NUM, size_t RHS_COL_NUM>
auto operator*(const Matrix<Type2, ROW_NUM, COL_NUM>& lhs, const Vector<Type, COL_NUM, RHS_COL_NUM>& rhs)
-> Vector<decltype(std::declval<Type2>() - std::declval<Type>()), ROW_NUM, RHS_COL_NUM>{

	return (lhs*static_cast<const Matrix<Type, COL_NUM, RHS_COL_NUM>&>(rhs));
}

template <typename Type, size_t size>
using ColVector = Vector< Type, size, 1 >;


typedef ColVector<double, 2> Vec2d;
typedef ColVector<int, 2> Vec2i;
typedef ColVector<double, 3> Vec3d;
typedef ColVector<int, 3> Vec3i;


#endif
