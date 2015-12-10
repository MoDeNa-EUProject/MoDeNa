
double determinant(double **in_matrix, int n);

double determinant(double **in_matrix, int n)
{
  int i, j, k;
  double **matrix;
  double det = 1;
  matrix = new double *[n];

  for (i = 0; i < n; i++)
  {
    matrix[i] = new double[n];
  }
  
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      matrix[i][j] = in_matrix[i][j];
    }
  }

  for (k = 0; k < n; k++)
  {
    if (matrix[k][k] == 0)
    {
      bool ok = false;

      for (j = k; j < n; j++)
      {
        if (matrix[j][k] != 0)
        {
          ok = true;
        }
      }

      if (!ok) {return 0;}

      for (i = k; i < n; i++)
      {
        std::swap (matrix[i][j], matrix[i][k]);
      }

      det = -det;
    }

    det *= matrix[k][k];

    if (k+1 < n)
    {
      for (i = k+1; i < n; i++)
      {
        for (j = k+1; j < n; j++)
        {
          matrix[i][j] = matrix[i][j] - matrix[i][k]*matrix[k][j]/matrix[k][k];
        }
      }
    }
  }

  // for (i = 0; i < n; i++)
  //   delete [] matrix[i];

  // delete [] matrix;

  return det;
}

