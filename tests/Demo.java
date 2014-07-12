package tests;

import utils.ArrayHelper;

public class Demo {

    private static final int ROUND_NUMBER = 3;

    public static void main(String[] args) {
        int n = 2;

        double[][] a = new double[n][n];
        double[][] b = new double[n][n];
        double[][] c = new double[n][n];

        // fillWithRandomStochasticMatrix(a);
        // fillWithRandomStochasticMatrix(b);

        a = new double[][] { { 1, 2 }, { 3, 4 } };
        b = new double[][] { { 0, 1 }, { 2, 0 } };

        printArray(a);
        System.out.println("*");
        printArray(b);
        System.out.println("=");

        matrixMult(a, b, c);

        printArray(c);
    }

    private static void matrixMult(double[][] a, double[][] b, double[][] c) {
        if (!dimensionsMatch(a, b, c))
            throw new IllegalArgumentException();

        for (int row = 0; row < a.length; ++row) {
            for (int col = 0; col < b[0].length; ++col) {
                c[row][col] = 0.0;
                for (int i = 0; i < a[0].length; ++i) {
                    c[row][col] += a[row][i] * b[i][col];
                }
            }
        }
    }

    private static boolean dimensionsMatch(double[][] a, double[][] b,
            double[][] c) {
        boolean result = isMatrix(a) && isMatrix(b) && isMatrix(c);
        if (result) {
            return a[0].length == b.length && c.length == a.length
                    && c[0].length == b[0].length;
        } else {
            return false;
        }
    }

    private static boolean isMatrix(double[][] m) {
        int leftCols = m[0].length;
        for (int i = 1; i < m.length; ++i) {
            if (m[i].length != leftCols) {
                return false;
            }
        }
        return true;
    }

    private static void printArray(double[][] m) {
        for (int i = 0; i < m.length; ++i) {
            System.out.println(ArrayHelper.toString(m[i], ROUND_NUMBER));
        }
    }

    private static void fillWithRandomStochasticMatrix(double[][] m) {
        for (int i = 0; i < m.length; ++i) {
            double valLeft = 1.0;
            for (int j = 0; j < m[i].length; ++j) {
                m[i][j] = Math.min(Math.random() * 2 / m[i].length, valLeft);
                valLeft -= m[i][j];
            }
            m[i][m[i].length - 1] += valLeft;
        }
    }

    private static void fillWithUniformDistribution(double[][] m) {
        for (int i = 0; i < m.length; ++i) {
            for (int j = 0; j < m[i].length; ++j) {
                m[i][j] = 1.0 / m[i].length;
            }
        }
    }
}
