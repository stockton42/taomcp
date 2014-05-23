package matrices;

import java.util.concurrent.RecursiveTask;

public abstract class Matrix {

    // determines to which size the strassen algorithm should use the standard
    // multiplication
    private static final int STRASSEN_MULT_LIMIT = 32;

    public abstract Matrix clone();

    public boolean multPossible(Matrix matrix) {
        return this.getCols() == matrix.getRows();
    }

    public abstract Matrix multWith(Matrix matrix);

    public abstract Matrix prlMultWith(Matrix matrix);

    public abstract double get(int row, int col);

    public abstract void put(double val, int row, int col);

    public abstract void del(int row, int col);

    public abstract int getRows();

    public abstract int getCols();

    public abstract Matrix getPart(int row1, int col1, int row2, int col2);

    public abstract Matrix getNewInstance(int rows, int cols);

    public abstract void add(Matrix mat);

    public Matrix cloneAdd(Matrix mat) {
        Matrix result = this.clone();
        result.add(mat);
        return result;
    }

    public abstract void sub(Matrix mat);

    public Matrix cloneSub(Matrix mat) {
        Matrix result = this.clone();
        result.sub(mat);
        return result;
    }

    public boolean hasSameDimensions(Matrix mat) {
        return this.getCols() == mat.getCols()
                && this.getRows() == mat.getRows();
    }

    public abstract void strassenMultThisWithInto(Matrix matrix, Matrix result);

    protected void strassenMultThisWithInto(Matrix matrix, Matrix result,
            boolean writeByRow) {
        // calculate size for strassen algorithm
        int size = getStrassenCalcSize(matrix);

        if (size < 1) {
            throw new IllegalStateException();
        } else if (size == 1) {
            result.put(this.get(0, 0) * matrix.get(0, 0), 0, 0);
        } else if (size <= STRASSEN_MULT_LIMIT) {
            multThisWithInto(matrix, result, writeByRow);
        } else if (size == 2) {
            double i, ii, iii, iv, v, vi, vii;
            double a11, a12, a21, a22, b11, b12, b21, b22;
            a11 = this.get(0, 0);
            a12 = this.get(0, 1);
            a21 = this.get(1, 0);
            a22 = this.get(1, 1);
            b11 = matrix.get(0, 0);
            b12 = matrix.get(0, 1);
            b21 = matrix.get(1, 0);
            b22 = matrix.get(1, 1);

            i = (a12 - a22) * (b21 + b22);
            ii = (a11 + a22) * (b11 + b22);
            iii = (a11 - a21) * (b11 + b12);
            iv = (a11 + a12) * b22;
            v = a11 * (b12 - b22);
            vi = a22 * (b21 - b11);
            vii = (a21 + a22) * b11;

            result.put(i + ii - iv + vi, 0, 0);
            result.put(iv + v, 0, 1);
            result.put(vi + vii, 1, 0);
            result.put(ii - iii + v - vii, 1, 1);

        } else {
            int outerHalfSize = size / 2;
            int innerHalfSize = outerHalfSize - 1;

            Matrix A11 = this.getPart(0, 0, innerHalfSize, innerHalfSize);
            Matrix A21 = this
                    .getPart(outerHalfSize, 0, size - 1, innerHalfSize);
            Matrix A12 = this
                    .getPart(0, outerHalfSize, innerHalfSize, size - 1);
            Matrix A22 = this.getPart(outerHalfSize, outerHalfSize, size - 1,
                    size - 1);

            Matrix B11 = matrix.getPart(0, 0, innerHalfSize, innerHalfSize);
            Matrix B21 = matrix.getPart(outerHalfSize, 0, size - 1,
                    innerHalfSize);
            Matrix B12 = matrix.getPart(0, outerHalfSize, innerHalfSize,
                    size - 1);
            Matrix B22 = matrix.getPart(outerHalfSize, outerHalfSize, size - 1,
                    size - 1);

            Matrix I = result.getNewInstance(outerHalfSize, outerHalfSize);
            (A12.cloneSub(A22)).strassenMultThisWithInto(B21.cloneAdd(B22), I,
                    writeByRow);
            Matrix II = result.getNewInstance(outerHalfSize, outerHalfSize);
            (A11.cloneAdd(A22)).strassenMultThisWithInto(B11.cloneAdd(B22), II,
                    writeByRow);
            Matrix III = result.getNewInstance(outerHalfSize, outerHalfSize);
            (A11.cloneSub(A21)).strassenMultThisWithInto(B11.cloneAdd(B12),
                    III, writeByRow);
            Matrix IV = result.getNewInstance(outerHalfSize, outerHalfSize);
            (A11.cloneAdd(A12)).strassenMultThisWithInto(B22, IV, writeByRow);
            Matrix V = result.getNewInstance(outerHalfSize, outerHalfSize);
            A11.strassenMultThisWithInto(B12.cloneSub(B22), V, writeByRow);
            Matrix VI = result.getNewInstance(outerHalfSize, outerHalfSize);
            A22.strassenMultThisWithInto(B21.cloneSub(B11), VI, writeByRow);
            Matrix VII = result.getNewInstance(outerHalfSize, outerHalfSize);
            (A21.cloneAdd(A22)).strassenMultThisWithInto(B11, VII, writeByRow);

            Matrix C11 = I;
            I.add(II);
            I.sub(IV);
            I.add(VI);

            Matrix C12 = IV;
            IV.add(V);

            Matrix C21 = VI;
            VI.add(VII);

            Matrix C22 = II;
            II.sub(III);
            II.add(V);
            II.sub(VII);

            result.pool(C11, C12, C21, C22);
        }
    }

    public abstract void pool(Matrix upLeft, Matrix upRight, Matrix downLeft,
            Matrix downRight);

    public boolean poolPossible(Matrix upLeft, Matrix upRight, Matrix downLeft,
            Matrix downRight) {
        return upLeft.hasSameDimensions(upRight)
                && upRight.hasSameDimensions(downLeft)
                && downLeft.hasSameDimensions(downRight)
                && upLeft.getRows() == upLeft.getCols();
    }

    private int getStrassenCalcSize(Matrix matrix) {
        int max = Math.max(
                this.getCols(),
                Math.max(this.getRows(),
                        Math.max(matrix.getCols(), matrix.getRows())));
        int newSize = 1;
        while (newSize < max) {
            newSize *= 2;
        }
        return newSize;
    }

    protected void prlStrassenMultThisWithInto(Matrix matrix, Matrix result) {
        // TODO
    }

    public abstract void winogradMultThisWithInto(Matrix matrix, Matrix result);

    protected void winogradMultThisWithInto(Matrix matrix, Matrix result,
            boolean writeByRow) {
        double[] A = new double[this.getRows()];
        double[] B = new double[matrix.getCols()];

        for (int leftRow = 0; leftRow < this.getRows(); ++leftRow) {
            A[leftRow] = 0;
            for (int leftCol = 0; leftCol < this.getCols() / 2; ++leftCol) {
                A[leftRow] += this.get(leftRow, 2 * leftCol + 1)
                        * this.get(leftRow, 2 * leftCol);
            }
        }
        for (int rightCol = 0; rightCol < matrix.getCols(); ++rightCol) {
            B[rightCol] = 0;
            for (int rightRow = 0; rightRow < matrix.getRows() / 2; ++rightRow) {
                B[rightCol] += matrix.get(2 * rightRow + 1, rightCol)
                        * matrix.get(2 * rightRow, rightCol);
            }
        }

        if (writeByRow) {
            for (int leftRow = 0; leftRow < this.getRows(); ++leftRow) {
                for (int rightCol = 0; rightCol < matrix.getCols(); ++rightCol) {
                    writeEntryWithTo(matrix, result, A, B, leftRow, rightCol);
                }
            }
        } else {
            for (int rightCol = 0; rightCol < matrix.getCols(); ++rightCol) {
                for (int leftRow = 0; leftRow < this.getRows(); ++leftRow) {
                    writeEntryWithTo(matrix, result, A, B, leftRow, rightCol);
                }
            }
        }
    }

    private void writeEntryWithTo(Matrix matrix, Matrix result, double[] A,
            double[] B, int leftRow, int rightCol) {
        double val = 0;
        for (int leftCol = 0; leftCol <= this.getCols() / 2; ++leftCol) {
            val += (this.get(leftRow, 2 * leftCol) + matrix.get(
                    2 * leftCol + 1, rightCol))
                    * (this.get(leftRow, 2 * leftCol + 1) + matrix.get(
                            2 * leftCol, rightCol));
        }
        result.put(val - A[leftRow] - B[rightCol], leftRow, rightCol);
    }

    public abstract double getMaxNorm();

    public abstract double get2Norm();

    public boolean equals(Matrix mat) {
        if (this.getCols() != mat.getCols() || this.getRows() != mat.getRows()) {
            return false;
        }

        for (int c = 0; c < this.getCols(); ++c) {
            for (int r = 0; r < this.getRows(); ++r) {
                if (this.get(r, c) != mat.get(r, c)) {
                    return false;
                }
            }
        }

        return true;
    }

    protected void multThisWithInto(Matrix matrix, Matrix result,
            boolean writeByRow) {
        if (writeByRow) {
            for (int thisRow = 0; thisRow < this.getRows(); ++thisRow) {
                for (int thatCol = 0; thatCol < matrix.getCols(); ++thatCol) {
                    writeEntryFromInto(matrix, result, thisRow, thatCol);
                }
            }
        } else {
            for (int thatCol = 0; thatCol < matrix.getCols(); ++thatCol) {
                for (int thisRow = 0; thisRow < this.getRows(); ++thisRow) {
                    writeEntryFromInto(matrix, result, thisRow, thatCol);
                }
            }
        }
    }

    private void writeEntryFromInto(Matrix matrix, Matrix result, int thisRow,
            int thatCol) {
        double entry = 0;
        for (int thisCol = 0; thisCol < this.getCols(); ++thisCol) {
            entry += this.get(thisRow, thisCol) * matrix.get(thisCol, thatCol);
        }
        result.put(entry, thisRow, thatCol);
    }

    protected void prlMultThisWithInto(Matrix matrix, Matrix[] result,
            int threads, boolean computeColumns) {
        ColumnProcessor[] workers = new ColumnProcessor[threads];

        for (int a = 0; a < threads; a++) {
            workers[a] = new ColumnProcessor(a, threads, this, matrix,
                    result[a], computeColumns);
            workers[a].fork();
        }
        for (int a = threads - 1; a >= 0; a--) {
            workers[a].join();
        }
    }

    public String toString() {
        StringBuilder result = new StringBuilder("");
        for (int i = 0; i < this.getRows(); i++) {
            for (int j = 0; j < this.getCols(); j++) {
                double entry = this.get(i, j);
                result.append(entry + "\t");
            }

            result.append("\n");
        }

        if (result.length() > 0) {
            // no last row
            result.deleteCharAt(result.length() - 1);
        }

        return result.toString();
    }

    public boolean isValidEntryLocation(int row, int col) {
        return row >= 0 && col >= 0 && row < this.getRows()
                && col < this.getCols();
    }

    private class ColumnProcessor extends RecursiveTask<Double> {
        private static final long serialVersionUID = 4831852382173221932L;
        private int id, threads;
        private Matrix left, right, target;
        private boolean computeColumns = false;

        private ColumnProcessor(int id, int threads, Matrix left, Matrix right,
                Matrix result, boolean computeColumns) {
            this.computeColumns = computeColumns;
            this.id = id;
            this.threads = threads;
            this.left = left;
            this.right = right;
            this.target = result;
        }

        @Override
        protected Double compute() {
            if (computeColumns) {
                // computes every (id + i*threads) column for i = 0, 1*threads,
                // 2*threads...
                for (int leftRow = 0; leftRow < left.getRows(); leftRow++) {
                    for (int rightCol = id; rightCol < right.getCols(); rightCol += threads) {
                        writeEntryAt(leftRow, rightCol);
                    }
                }
            } else {
                // computes every (id + i*threads) row for i = 0, 1*threads,
                // 2*threads...
                for (int leftRow = id; leftRow < left.getRows(); leftRow += threads) {
                    for (int rightCol = 0; rightCol < right.getCols(); rightCol++) {
                        writeEntryAt(leftRow, rightCol);
                    }
                }
            }

            // no result needed
            return 0.0;
        }

        private void writeEntryAt(int leftRow, int rightCol) {
            double temp = 0;
            for (int c = 0; c < left.getCols(); c++) {
                temp += left.get(leftRow, c) * right.get(c, rightCol);
            }
            target.put(temp, leftRow, rightCol);
        }
    }
}
