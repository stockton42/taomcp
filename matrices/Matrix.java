package matrices;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.concurrent.RecursiveTask;

public abstract class Matrix {

    public static final boolean SHOW_MATRICES_WHEN_NOT_EQUAL = true;

    // set up number format for toString() method
    private static final DecimalFormatSymbols ds = new DecimalFormatSymbols();
    private static final DecimalFormat df;
    static {
        df = new DecimalFormat("#0.00", ds);
    }

    /**
     * Determines to which size the Strassen algorithm should use the standard
     * multiplication.
     */
    private static final int STRASSEN_MATRIX_SIZE_LIMIT = 32;

    /**
     * The Strassen algorithm will calculate the submatrices in separate threads
     * until the recursion reaches this level.
     */
    private static final int STRASSEN_PARALLEL_LEVEL_LIMIT = 1;

    private static final boolean SHOW_MODIFICATIONS = true;

    public abstract Matrix clone();

    public boolean multPossible(Matrix matrix) {
        return this.getCols() == matrix.getRows();
    }

    protected abstract Matrix multWith(Matrix matrix);

    protected abstract Matrix prlMultWith(Matrix matrix);

    public abstract double get(int row, int col);

    public abstract void put(double val, int row, int col);

    public abstract void del(int row, int col);

    public abstract int getRows();

    public abstract int getCols();

    public abstract Matrix getPart(int row1, int col1, int row2, int col2);

    public abstract Matrix getNewInstance(int rows, int cols);

    /**
     * Checks if this matrix is nonnegative.
     * 
     * @return is true if every entry of this matrix is greater than or equal to
     *         zero.
     */
    public abstract boolean isNonNegative();

    public void setNegativeEntriesToZero() {
        setNegativeEntriesToZero(SHOW_MODIFICATIONS);
    }

    protected abstract void setNegativeEntriesToZero(boolean showModifications);

    public abstract void add(Matrix mat);

    public Matrix getZero() {
        return getNewInstance(this.getRows(), this.getCols());
    }

    public Matrix getOne() {
        Matrix result = getNewInstance(this.getRows(), this.getCols());

        for (int i = 0; i < Math.min(this.getRows(), this.getCols()); ++i) {
            result.put(1, i, i);
        }

        return result;
    }

    public abstract void stabilizeRowsTo(double stabilizeRowsTo);

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

    /**
     * Multiplies a matrix with another matrix using a given multiplication
     * algorithm. The result will be returned, this matrix will NOT be modified.
     * 
     * @param matrix
     * @param multType
     * @return
     */
    public Matrix multWith(Matrix matrix, MatrixMultType multType) {
        switch (multType) {
        case NAIVE:
            return this.multWith(matrix);
        case PARALLEL_NAIVE:
            return this.prlMultWith(matrix);
        case PARALLEL_STRASSEN_NAIVE_HYBRID:
            return this.prlStrassenMultThisWith(matrix);
        case STRASSEN_NAIVE_HYBRID:
            return this.strassenMultThisWith(matrix);
        case WINOGRAD:
            return this.winogradMultThisWith(matrix);
        default:
            throw new IllegalArgumentException();
        }
    }

    protected abstract Matrix strassenMultThisWith(Matrix matrix);

    protected void strassenMultThisWithInto(Matrix matrix, Matrix result,
            boolean writeByRow) {
        strassenMultThisWithInto(matrix, result, writeByRow,
                STRASSEN_MATRIX_SIZE_LIMIT);
    }

    private void strassenMultThisWithInto(Matrix matrix, Matrix result,
            boolean writeByRow, final int level) {
        // calculate size for strassen algorithm
        int size = getStrassenCalcSize(matrix);

        if (size < 1) {
            throw new IllegalStateException();
        } else if (size == 1) {
            result.put(this.get(0, 0) * matrix.get(0, 0), 0, 0);
        } else if (size <= STRASSEN_MATRIX_SIZE_LIMIT) {
            multThisWithInto(matrix, result, writeByRow);
        } else if (size == 2) {
            strassenMultForTwoTimesTwo(matrix, result);
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
            Matrix II = result.getNewInstance(outerHalfSize, outerHalfSize);
            Matrix III = result.getNewInstance(outerHalfSize, outerHalfSize);
            Matrix IV = result.getNewInstance(outerHalfSize, outerHalfSize);
            Matrix V = result.getNewInstance(outerHalfSize, outerHalfSize);
            Matrix VI = result.getNewInstance(outerHalfSize, outerHalfSize);
            Matrix VII = result.getNewInstance(outerHalfSize, outerHalfSize);

            RecursiveTask<Double> recTaskForI = null;
            RecursiveTask<Double> recTaskForII = null;
            RecursiveTask<Double> recTaskForIII = null;
            RecursiveTask<Double> recTaskForIV = null;
            RecursiveTask<Double> recTaskForV = null;
            RecursiveTask<Double> recTaskForVI = null;

            if (level < STRASSEN_PARALLEL_LEVEL_LIMIT) {
                recTaskForI = new RecursiveTask<Double>() {
                    private static final long serialVersionUID = -9073807035893659703L;

                    @Override
                    protected Double compute() {
                        (A12.cloneSub(A22)).strassenMultThisWithInto(
                                B21.cloneAdd(B22), I, writeByRow, level + 1);
                        return 0.0; // no result needed
                    }
                };
                recTaskForI.fork();

                recTaskForII = new RecursiveTask<Double>() {
                    private static final long serialVersionUID = 6290324240156453866L;

                    @Override
                    protected Double compute() {
                        (A11.cloneAdd(A22)).strassenMultThisWithInto(
                                B11.cloneAdd(B22), II, writeByRow, level + 1);
                        return 0.0; // no result needed
                    }
                };
                recTaskForII.fork();

                recTaskForIII = new RecursiveTask<Double>() {
                    private static final long serialVersionUID = -6281873677476022429L;

                    @Override
                    protected Double compute() {
                        (A11.cloneSub(A21)).strassenMultThisWithInto(
                                B11.cloneAdd(B12), III, writeByRow, level + 1);
                        return 0.0; // no result needed
                    }
                };
                recTaskForIII.fork();

                recTaskForIV = new RecursiveTask<Double>() {
                    private static final long serialVersionUID = -1839210672357109597L;

                    @Override
                    protected Double compute() {
                        (A11.cloneAdd(A12)).strassenMultThisWithInto(B22, IV,
                                writeByRow, level + 1);
                        return 0.0; // no result needed
                    }
                };
                recTaskForIV.fork();

                recTaskForV = new RecursiveTask<Double>() {
                    private static final long serialVersionUID = -4948362458881227395L;

                    @Override
                    protected Double compute() {
                        A11.strassenMultThisWithInto(B12.cloneSub(B22), V,
                                writeByRow, level + 1);
                        return 0.0; // no result needed
                    }
                };
                recTaskForV.fork();

                recTaskForVI = new RecursiveTask<Double>() {
                    private static final long serialVersionUID = 7273900793962698628L;

                    @Override
                    protected Double compute() {
                        A22.strassenMultThisWithInto(B21.cloneSub(B11), VI,
                                writeByRow, level + 1);
                        return 0.0; // no result needed
                    }
                };
                recTaskForVI.fork();
            } else {
                (A12.cloneSub(A22)).strassenMultThisWithInto(B21.cloneAdd(B22),
                        I, writeByRow, level + 1);
                (A11.cloneAdd(A22)).strassenMultThisWithInto(B11.cloneAdd(B22),
                        II, writeByRow, level + 1);
                (A11.cloneSub(A21)).strassenMultThisWithInto(B11.cloneAdd(B12),
                        III, writeByRow, level + 1);
                (A11.cloneAdd(A12)).strassenMultThisWithInto(B22, IV,
                        writeByRow, level + 1);
                A11.strassenMultThisWithInto(B12.cloneSub(B22), V, writeByRow,
                        level + 1);
                A22.strassenMultThisWithInto(B21.cloneSub(B11), VI, writeByRow,
                        level + 1);
            }

            // do the last calculation in this thread
            (A21.cloneAdd(A22)).strassenMultThisWithInto(B11, VII, writeByRow,
                    STRASSEN_PARALLEL_LEVEL_LIMIT);

            if (level < STRASSEN_PARALLEL_LEVEL_LIMIT) {
                recTaskForVI.join();
                recTaskForV.join();
                recTaskForIV.join();
                recTaskForIII.join();
                recTaskForII.join();
                recTaskForI.join();
            }

            collectSubMatrices(result, I, II, III, IV, V, VI, VII);
        }
    }

    private void collectSubMatrices(Matrix result, Matrix I, Matrix II,
            Matrix III, Matrix IV, Matrix V, Matrix VI, Matrix VII) {
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

    private void strassenMultForTwoTimesTwo(Matrix matrix, Matrix result) {
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
    }

    protected abstract void pool(Matrix upLeft, Matrix upRight,
            Matrix downLeft, Matrix downRight);

    protected boolean poolPossible(Matrix upLeft, Matrix upRight,
            Matrix downLeft, Matrix downRight) {
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

    protected abstract Matrix prlStrassenMultThisWith(Matrix matrix);

    protected void prlStrassenMultThisWithInto(Matrix matrix, Matrix result,
            boolean writeByRow) {
        strassenMultThisWithInto(matrix, result, writeByRow, 0);
    }

    protected abstract Matrix winogradMultThisWith(Matrix matrix);

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

    public abstract double getNorm(MatrixNorm norm);

    public boolean equals(Matrix mat) {
        if (this.getCols() != mat.getCols() || this.getRows() != mat.getRows()) {
            return false;
        }

        for (int c = 0; c < this.getCols(); ++c) {
            for (int r = 0; r < this.getRows(); ++r) {
                if (this.get(r, c) != mat.get(r, c)) {

                    if (SHOW_MATRICES_WHEN_NOT_EQUAL) {
                        // System.err.println("\n" + this);
                        System.err.println("\nERROR AT ENTRY (" + r + ", " + c
                                + "): " + this.get(r, c) + " VS "
                                + mat.get(r, c) + ":\t"
                                + (this.get(r, c) - mat.get(r, c)) + "\n");
                        // System.err.println(mat);
                    }

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
            int threads, boolean writeByRow) {
        ColumnProcessor[] workers = new ColumnProcessor[threads];

        for (int a = 0; a < threads; a++) {
            workers[a] = new ColumnProcessor(a, threads, this, matrix,
                    result[a], writeByRow);
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
                result.append(df.format(entry) + "\t");
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
        private boolean writeByRow = false;

        private ColumnProcessor(int id, int threads, Matrix left, Matrix right,
                Matrix result, boolean writeByRow) {
            this.writeByRow = writeByRow;
            this.id = id;
            this.threads = threads;
            this.left = left;
            this.right = right;
            this.target = result;
        }

        @Override
        protected Double compute() {
            if (!writeByRow) {
                // computes every (id + i*threads) column for i = 0, 1*threads,
                // 2*threads...
                for (int rightCol = id; rightCol < right.getCols(); rightCol += threads) {
                    for (int leftRow = 0; leftRow < left.getRows(); leftRow++) {
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

            return 0.0; // no result needed
        }

        private void writeEntryAt(int leftRow, int rightCol) {
            double temp = 0;
            for (int c = 0; c < left.getCols(); c++) {
                temp += left.get(leftRow, c) * right.get(c, rightCol);
            }
            target.put(temp, leftRow, rightCol);
        }
    }

    /**
     * Returns the minimal positive entry of this matrix. If all entries are
     * equal to zero, zero is returned. The behavior is well-defined if and only
     * if this matrix is a stochastic matrix.
     * 
     * @return a positive double value if there is at least one positive entry,
     *         zero otherwise.
     */
    public abstract double getMinimalPositiveEntry();
}
