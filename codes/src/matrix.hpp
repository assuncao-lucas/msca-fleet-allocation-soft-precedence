#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>
#include <cassert>

template <class T>
class Matrix final
{
public:
	explicit Matrix() = default;
	explicit Matrix(const int, const int);
	explicit Matrix(const int, const int, T);

	Matrix(const Matrix<T> &src);
	Matrix(const Matrix<T> &&src);

	virtual ~Matrix();

	Matrix<T> *Transposed();
	Matrix<T> operator+(Matrix<T> &) const;
	T *operator[](int) const;
	// copy assignment operator.
	Matrix<T> &operator=(const Matrix<T> &src);
	// move assignment operator.
	Matrix<T> &operator=(Matrix<T> &&src);
	int lines() const { return lines_; }
	int columns() const { return columns_; }
	T **matrix() const { return matrix_; }
	void set_matrix(T **matrix) { matrix_ = matrix; }
	template <class T2>
	friend std::ostream &operator<<(std::ostream &out, const Matrix<T2> &matrix);

protected:
	int lines_ = 0;
	int columns_ = 0;
	T **matrix_ = nullptr;
	void Init(const int lines, const int columns);
	void Init(const int lines, const int columns, T);
	void Destroy();
	friend void swap(Matrix<T> &first, Matrix<T> &second) noexcept
	{
		for (int i = 1; i < lines; ++i)
			std::swap(first.matrix_[i], second.matrix_[i]);
	}
};

template <class T>
void Matrix<T>::Init(const int lines, const int columns)
{
	lines_ = lines;
	columns_ = columns;

	matrix_ = new T *[lines];
	(matrix_)[0] = new T[lines * columns];
	for (int i = 1; i < lines; ++i)
		(matrix_)[i] = (matrix_)[0] + i * columns;
}

template <class T>
void Matrix<T>::Init(const int lines, const int columns, T initial_value)
{
	Init(lines, columns);

	for (int i = 0; i < lines; ++i)
		for (int j = 0; j < columns; ++j)
			(*this)[i][j] = initial_value;
}

template <typename T>
Matrix<T>::Matrix(const int lines, const int columns)
{
	Init(lines, columns);
}

template <typename T>
Matrix<T>::Matrix(const int lines, const int columns, T initial_value)
{
	Init(lines, columns, initial_value);
}

template <typename T>
Matrix<T>::Matrix(const Matrix<T> &src)
{
	Init(src.lines_, src.columns_);

	for (int i = 0; i < src.lines_; ++i)
		for (int j = 0; j < src.columns_; ++j)
			(*this)[i][j] = src[i][j];
}

template <typename T>
Matrix<T>::Matrix(const Matrix<T> &&src)
{
	assert(lines_ == src.lines_ && columns_ == src.columns_);
	swap(*this, src);
}

template <typename T>
Matrix<T> *Matrix<T>::Transposed()
{
	Matrix<T> *transposed = new Matrix<T>(columns_, lines_);
	for (int i = 0; i < lines_; ++i)
		for (int j = 0; j < columns_; ++j)
			(*transposed)[j][i] = (*this)[i][j];

	return transposed;
}

template <typename T>
Matrix<T>::~Matrix()
{
	Destroy();
}

template <typename T>
void Matrix<T>::Destroy()
{
	if (matrix_)
	{
		delete[] matrix_[0];
		delete[] matrix_;
		matrix_ = nullptr;
	}
}

template <class T>
Matrix<T> Matrix<T>::operator+(Matrix<T> &m) const
{
	if ((lines_ != m.lines_) || (columns_ != m.columns_))
		throw 1;

	Matrix<T> sum_matrix(lines_, columns_);

	for (int i = 0; i < lines_; ++i)
		for (int j = 0; j < columns_; ++j)
			sum_matrix[i][j] = (*this)[i][j] + m[i][j];

	return sum_matrix;
}

template <class T>
Matrix<T> &Matrix<T>::operator=(Matrix<T> &&other) // move assignment
{
	if (this != &other)
	{
		Destroy(); // destroy storage in this.
		matrix_ = other.matrix();
		other.set_matrix(nullptr);
	}

	return *this;
}

template <class T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &other) // copy assignment
{
	if (this != &other)
	{ // self-assignment check expected
		if ((lines_ != other.lines()) || (columns_ != other.columns()))
		{
			Destroy();							  // destroys storage in this
			Init(other.lines(), other.columns()); // creates storage in this
		}
		/* copies data from other's storage to this storage */

		for (int i = 0; i < other.lines(); ++i)
		{
			for (int j = 0; j < other.columns(); ++j)
			{
				matrix_[i][j] = other[i][j];
			}
		}
	}
	return *this;
}

template <class T>
T *Matrix<T>::operator[](int i) const
{
	if (i >= lines_)
	{
		std::cout << "line " << i << " is greater than dimension " << lines_ << std::endl;
		getchar();
		getchar();
	}
	return (this->matrix_)[i];
}

template <class T>
std::ostream &operator<<(std::ostream &out, Matrix<T> &matrix)
{
	for (int i = 0; i < matrix.lines(); ++i)
	{
		for (int j = 0; j < matrix.columns(); ++j)
		{
			out << matrix[i][j] << " ";
		}
		out << std::endl;
	}
	return out;
}

#endif