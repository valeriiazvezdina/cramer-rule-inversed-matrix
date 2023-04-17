#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

const int dimension = 4; // the dimension of the matrix

void print_matrix(const vector<vector<double>>& m) // function for printing the matrix
{
    for (const auto &i : m)
    {
        cout << "|"; // vizualizing like in a determinant braces

        for (const auto &j : i)
        {
            cout << " " << j << " ";
        }

        cout << "|" << endl;
    }
}

vector<vector<double>> initMatrix() // function for initialization the matrix which asks values from user
{
    vector<vector<double>> m (dimension, vector<double>(dimension));

    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            cout << "enter the "<< j + 1 << " element of the " << i + 1 << " row of the matrix" << endl;
            cin >> m[i][j];
        }
    }

    // det = 117
//    m.push_back({ 6, 7, 1 });
//    m.push_back({ 3, 4, 2 });
//    m.push_back({ 9, -2, 3 });

    // det = 6
//    m.push_back({ 1, 2, 3 });
//    m.push_back({ 1, 0, 1 });
//    m.push_back({ 4, 5, 6 });

//    m.push_back({ 1, -1, 0 });
//    m.push_back({ 0, 2, 0 });
//    m.push_back({ 0, 1, -1 });

    return m;
}

vector<double> getRightSide() // function which asks the right side vector from the user
{
    vector<double> side (dimension);

    for (int i = 0; i < dimension; i++)
    {
        cout << "enter the "<< i + 1 << " element of the right side vector" << endl;
        cin >> side[i];
    }

    return side;
}

double calculate_det(vector<vector<double>>& matrix, int dim) // function for calculating the determinant
{
    // some default values
    vector<vector<double>> m(matrix);
    double det = 1;
    double delta = 0;
    int sgn = 1;

    if (dim == 1) return m[0][0]; // if dimension is 1 there is nothing to calculate

    if (dim == 2) return (m[0][0] * m[1][1] - m[1][0] * m[0][1]);

    for (int i = 0; i < dim; i++)
    {
        if (m[i][i] == 0)
        {
            // checking if the first element of the next row is 0 or no
            for (int k = i + 1; k < dim; k++)
            {
                if (k >= dim) return 0; // if the k is out of dimension, so that matrix is singular

                if (m[k][i] != 0) // if yes, then skipping this row and check next after current one (k-th), if no, just swapping the rows
                {
                    // swapping the rows
                    auto temp = m[i];
                    m[i] = m[k];
                    m[k] = temp;
                    temp.clear();
                    sgn *= -1;
                    break;
                }
            }
        }

        det *= m[i][i]; // calculating the product of determinant

        for (int j = i + 1; j < dim; j++)
        {
            delta = m[j][i] / m[i][i]; // making row 0

            for (int k = i; k < dim; k++)
            {
                m[j][k] = m[j][k] - delta * m[i][k];
            }
        }
    }

    if (isnan(det) || det == 0.0) // escaping the errors of calculating by making them 0
    {
        return 0.0;
    }

    return det * sgn;
}

double cofactor(vector<vector<double>> m, int dim, int i, int j) // function for calculating the particular cofactor matrix and returns the value of its determinant
{
    vector<vector<double>> subm (dim - 1, vector<double>(dim - 1)); // sub matrix which is smaller by 1 row and 1 row than defualt matrix

    // some useful values
    double d = 1;
    int row = 0;
    int col = 0;

    // deleting the i row and the j column
    for (int ii = 0; ii < dim; ii++)
    {
        for (int jj = 0; jj < dim; jj++)
        {
            if (ii != i && jj != j)
            {
                subm[row][col++] = m[ii][jj];

                if (col == dim - 1)
                {
                    col = 0;
                    row++;
                }
            }
        }
    }

    double det = calculate_det(subm, dim - 1); // sub determinant

    if (isnan(det) || det == 0.0) // escaping the calculating errors
    {
        return 0.0;
    }
    else
    {
        d = (((i + j) % 2 == 0) ? 1 : -1) * det; // main calculating
    }

    return d;
}

vector<vector<double>> adjugated(vector<vector<double>> m, int dim)
{
    vector<vector<double>> adjM (dim, vector<double>(dim)); // function that makes en adjugated matrix using the cofactor and transpose this matrix

    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            adjM[j][i] = cofactor(m, dim, i, j); // initializing and transpose
        }
    }

    return adjM;
}

vector<vector<double>> inversed(vector<vector<double>> m, int dim) // function that inverse the matrix using agjugated matrix
{
    double det = calculate_det(m, dim);

    if (det == 0) // if determinant is 0 there is nothing to calculate
    {
        cout << "determinant is zero, so matrix is sungular and the inverse matrix does nor exist." << endl;
        exit(-1);
    }

    vector<vector<double>> adj = adjugated(m, dim);

    double invDet = 1 / det;

    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            adj[i][j] *= invDet;
        }
    }

    return adj;
}

void cramer_rule(vector<vector<double>> m, int dim, vector<double> right_side) // solving the equations using cramer's rule
{
    double detA = calculate_det(m, dim);

    if (detA == 0) // if determinant is 0 there is nothing to calculate
    {
        cout << "determinant is zero, cannot solve using Cramer's Rule." << endl;
        return;
    }

    for (int i = 0; i < dim; i++)
    {
        vector<vector<double>> temp(dim, vector<double>(dim));

        // making the temporary matrix like a copy of default matrix
        for (int j = 0; j < dim; j++)
        {
            for (int k = 0; k < dim; k++)
            {
                temp[j][k] += m[j][k];
            }
        }

        // changing the column to the column of the right side vector
        for (int j = 0; j < dim; j++)
        {
            temp[j][i] = right_side[j];
        }

        double detX = calculate_det(temp, dim); // determinant of the temporary matrix with changed column

        cout << "x" << i+1 << " = " << detX/detA << endl; // printing the result of calculating
    }
}

int main()
{
    cout << setprecision(2);

    vector<vector<double>> matrix = initMatrix();
    cout << "====== matrix ======" << endl;
    print_matrix(matrix);
    vector<vector<double>> inv = inversed(matrix, dimension);
    cout << "==== inverse matrix ====" << endl;
    print_matrix(inv);
    vector<double> right_side = getRightSide();
    cout << "==== result ====" << endl;
    cramer_rule(matrix, dimension, right_side);

    return 0;
}
