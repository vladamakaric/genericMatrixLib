#ifndef MATRIX_H
#define MATRIX_H

#include <array>
#include <type_traits>
#include <ostream>

template <typename Type, size_t ROW_NUM, size_t COL_NUM>
class Matrix{
protected:
	std::array<Type, ROW_NUM*COL_NUM> data;
	
	const Type& elem(int indx)const{
		return data[indx];
	}

	Type& elem(int indx){
		return const_cast<Type&>(static_cast<const Matrix<Type, ROW_NUM, COL_NUM>&>(*this).elem(indx));
	}

	const Type& elem(int row, int col)const{
		return data[row*COL_NUM + col];
	}

	Type& elem(int row, int col){
		return const_cast<Type&>(static_cast<const Matrix<Type, ROW_NUM, COL_NUM>&>(*this).elem(row, col));
	}
public:
	template <typename... T>
	Matrix(T... ts) : data({ ts... }) { }

	void setAll(Type val) {
		data.assign(val);
	}

	void setIdentity(){
		setAll(0);

		for (int i = 0; i < COL_NUM; i++)
			elem(i,i) = 1;
	}

	size_t getRowNum() const { return ROW_NUM; };
	size_t getColNum() const { return COL_NUM; };
	size_t getElemNum() const { return COL_NUM*ROW_NUM; };

	Type& operator()(int row, int col){
		return elem(row, col);
	}

	const Type& operator()(int row, int col) const{
		return elem(row, col);
	}

	const Type& operator()(int indx) const{
		return data[indx];
	}

	Type& operator()(int indx){
		return data[indx];
	}

	template <typename Type2>
	Matrix<Type, ROW_NUM, ROW_NUM>& operator*=(const Matrix<Type2, ROW_NUM, ROW_NUM>& rhs){
		Matrix temp(*this);
		MultiplyMatrices(*this, temp, rhs);
		return *this;
	}

	template <typename Type2>
	Matrix<Type, ROW_NUM, COL_NUM>& operator*=(const Type2& rhs){
		for (int i = 0; i < getElemNum(); i++)
			elem(i) *= rhs;

		return *this;
	}

	Matrix<Type, ROW_NUM, COL_NUM> operator-() const {

		Matrix<Type, ROW_NUM, COL_NUM> temp(*this);
		for (int i = 0; i < getElemNum(); i++)
			temp(i) *= -1;

		return temp;
	}


};

template <typename DestT, typename LhsT, typename RhsT,
	size_t DEST_ROW_NUM,
	size_t DEST_COL_NUM,
	size_t LHS_COL_NUM>
	static void MultiplyMatrices(Matrix<DestT, DEST_ROW_NUM, DEST_COL_NUM>& dest,
	const Matrix<LhsT, DEST_ROW_NUM, LHS_COL_NUM>& lhs,
	const Matrix<RhsT, LHS_COL_NUM, DEST_COL_NUM>& rhs) {

		for (int i = 0; i < DEST_ROW_NUM; i++){
			for (int j = 0; j < DEST_COL_NUM; j++){
				dest(i, j) = 0;

				for (int k = 0; k < LHS_COL_NUM; k++)
					dest(i, j) += lhs(i, k) * rhs(k, j);
			}
		}
}

template <typename Type1, typename Type2, size_t ROW_NUM, size_t COL_NUM>
Matrix<Type1, ROW_NUM, COL_NUM> operator/(const Matrix<Type1, ROW_NUM, COL_NUM>& lhs, const Type2& rhs) {

	Matrix<Type1, ROW_NUM, COL_NUM> temp(lhs);
	for (int i = 0; i < getElemNum(); i++)
		temp(i) /= rhs;

	return temp;
}


template <typename Type2, typename Type, size_t ROW_NUM, size_t COL_NUM, size_t RHS_COL_NUM>
auto operator*(const Matrix<Type2, ROW_NUM, COL_NUM>& lhs, const Matrix<Type, COL_NUM, RHS_COL_NUM>& rhs)
-> Matrix<decltype(std::declval<Type2>() - std::declval<Type>()), ROW_NUM, RHS_COL_NUM>{

	Matrix<decltype(std::declval<Type2>() - std::declval<Type>()), ROW_NUM, RHS_COL_NUM> temp;

	MultiplyMatrices(temp, lhs, rhs);

	return temp;
}
//template <typename Type2, typename Type, size_t ROW_NUM, size_t COL_NUM, size_t RHS_COL_NUM>
//auto operator*(const Matrix<Type2, ROW_NUM, COL_NUM>& lhs, const Matrix<Type, COL_NUM, RHS_COL_NUM>& rhs){
//
//	Matrix<decltype(std::declval<Type2>() - std::declval<Type>()), ROW_NUM, RHS_COL_NUM> temp;
//
//	MultiplyMatrices(temp, lhs, rhs);
//
//	return temp;
//}

template <typename MatType, typename DatType, size_t ROW_NUM, size_t COL_NUM>
Matrix<MatType, ROW_NUM, COL_NUM> operator*(const DatType& rhs, const Matrix<MatType, ROW_NUM, COL_NUM> &lhs){

	Matrix<MatType, ROW_NUM, COL_NUM> temp(lhs);
	for (int i = 0; i < temp.getElemNum(); i++)
		temp(i) *= rhs;

	return temp;
}

template <typename MatType, typename DatType, size_t ROW_NUM, size_t COL_NUM>
Matrix<MatType, ROW_NUM, COL_NUM> operator*(const Matrix<MatType, ROW_NUM, COL_NUM> &lhs, const DatType& rhs){

	Matrix<MatType, ROW_NUM, COL_NUM> temp(lhs);
	for (int i = 0; i < temp.getElemNum(); i++)
		temp(i) *= rhs;

	return temp;
}

template <typename Type1, typename Type2, size_t ROW_NUM, size_t COL_NUM> 
Matrix<decltype(std::declval<Type2>() - std::declval<Type1>()), ROW_NUM, COL_NUM> 
	operator+(const Matrix<Type1, ROW_NUM, COL_NUM> &rhs, const Matrix<Type2, ROW_NUM, COL_NUM> &lhs){

	Matrix<decltype(std::declval<Type2>() - std::declval<Type1>()), ROW_NUM, COL_NUM> temp;
	for (unsigned i = 0; i < temp.getElemNum(); i++)
		temp(i) = lhs(i) + rhs(i);

	return temp;
}

template <typename Type1, typename Type2, size_t ROW_NUM, size_t COL_NUM>
Matrix<decltype(std::declval<Type2>() - std::declval<Type1>()), ROW_NUM, COL_NUM>
operator-(const Matrix<Type1, ROW_NUM, COL_NUM> &rhs, const Matrix<Type2, ROW_NUM, COL_NUM> &lhs){

	return lhs + (-rhs);
}

template <typename Type, size_t ROW_NUM, size_t COL_NUM>
std::ostream &operator<< (std::ostream &out, const Matrix<Type, ROW_NUM, COL_NUM> &matrix){

	for (int i = 0; i < ROW_NUM; i++){
		out << "|";

		for (int j = 0; j < COL_NUM; j++)
			out << " " << matrix(i,j);

		out << " |" << std::endl;
	}

	return out;
}

typedef Matrix<double, 3, 3> Mat33d;
typedef Matrix<double, 2, 2> Mat22d;

#endif
