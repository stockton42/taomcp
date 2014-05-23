package matrices;

import java.util.Arrays;

public class ArrayMatrix extends Matrix {

    // calculate complete rows first, increases speed by 100 %
    private static final boolean WRITE_BY_ROW = true;

    private static final int DEFAULT_VALUE = 0;
    private static final int NUMBER_OF_THREADS = 8;

    private double[][] content;

    public ArrayMatrix(int rows, int cols) {
        this.content = new double[rows][cols];
    }

    public ArrayMatrix(double[][] content, boolean copy) {
        if (!isValidMatrix(content)) {
            throw new IllegalArgumentException();
        }

        if (copy) {
            this.content = new double[content.length][content[0].length];
            for (int i = 0; i < content.length; ++i) {
                this.content[i] = Arrays.copyOf(content[i], content[0].length);
            }
        } else {
            this.content = content;
        }
    }

    private boolean isValidMatrix(double[][] content) {
        if (content == null)
            return false;
        int rows = content.length;

        if (rows == 0 || content[0] == null) {
            return false;
        }

        int cols = content[0].length;

        for (int i = 1; i < content.length; ++i) {
            if (cols != content[i].length) {
                return false;
            }
        }

        return true;
    }

    @Override
    public double get(int row, int col) {
        if (!isValidEntryLocation(row, col))
            return DEFAULT_VALUE;

        return content[row][col];
    }

    @Override
    public void put(double val, int row, int col) {
        if (!isValidEntryLocation(row, col))
            return;

        content[row][col] = val;
    }

    @Override
    public void del(int row, int col) {
        if (!isValidEntryLocation(row, col))
            return;

        content[row][col] = DEFAULT_VALUE;
    }

    @Override
    public int getRows() {
        return content.length;
    }

    @Override
    public int getCols() {
        return content[0].length;
    }

    @Override
    public Matrix multWith(Matrix matrix) {
        if (!multPossible(matrix)) {
            throw new IllegalArgumentException();
        }

        ArrayMatrix result = new ArrayMatrix(this.getRows(), matrix.getCols());

        multThisWithInto(matrix, result, WRITE_BY_ROW);

        return result;
    }

    @Override
    public Matrix prlMultWith(Matrix matrix) {
        if (!multPossible(matrix)) {
            throw new IllegalArgumentException();
        }

        ArrayMatrix result = new ArrayMatrix(this.getRows(), matrix.getCols());
        ArrayMatrix[] matrices = new ArrayMatrix[NUMBER_OF_THREADS];
        for (int i = 0; i < NUMBER_OF_THREADS; ++i) {
            matrices[i] = result;
        }

        prlMultThisWithInto(matrix, matrices, NUMBER_OF_THREADS, WRITE_BY_ROW);

        return result;
    }

    @Override
    public Matrix getNewInstance(int rows, int cols) {
        return new ArrayMatrix(rows, cols);
    }

    @Override
    public Matrix getPart(int row1, int col1, int row2, int col2) {
        if (row1 > row2) {
            int temp = row2;
            row2 = row1;
            row1 = temp;
        }
        if (col1 > col2) {
            int temp = col1;
            col1 = col2;
            col2 = temp;
        }

        ArrayMatrix result = new ArrayMatrix(row2 - row1 + 1, col2 - col1 + 1);

        if (row1 < 0 || col1 < 0 || row2 < 0 || col2 < 0
                || col1 >= this.getCols() || row1 >= this.getRows()) {
            // invalid arguments: return empty matrix
            return result;
        }

        row2 = Math.min(this.getRows() - 1, row2);
        col2 = Math.min(this.getCols() - 1, col2);

        for (int row = row1; row <= row2; ++row) {
            for (int col = col1; col <= col2; ++col) {
                result.content[row - row1][col - col1] = this.content[row][col];
            }
        }

        return result;
    }

    @Override
    public void add(Matrix mat) {
        addSub(mat, true);
    }

    @Override
    public void sub(Matrix mat) {
        addSub(mat, false);
    }

    private void addSub(Matrix mat, boolean add) {
        if (hasSameDimensions(mat)) {
            for (int row = 0; row < this.getRows(); ++row) {
                for (int col = 0; col < this.getCols(); ++col) {
                    if (add) {
                        content[row][col] += mat.get(row, col);
                    } else {
                        content[row][col] -= mat.get(row, col);
                    }
                }
            }
        } else {
            throw new IllegalArgumentException();
        }
    }

    @Override
    public Matrix clone() {
        return new ArrayMatrix(this.content, true);
    }

    @Override
    public void pool(Matrix upLeft, Matrix upRight, Matrix downLeft,
            Matrix downRight) {
        if (!poolPossible(upLeft, upRight, downLeft, downRight)) {
            throw new IllegalArgumentException();
        }
        ArrayMatrix c11 = (ArrayMatrix) upLeft;
        ArrayMatrix c12 = (ArrayMatrix) upRight;
        ArrayMatrix c21 = (ArrayMatrix) downLeft;
        ArrayMatrix c22 = (ArrayMatrix) downRight;

        int halfSize = upLeft.getRows();
        ArrayMatrix leftCopyMat = c11;
        ArrayMatrix rightCopyMat = c12;
        int rowShift = 0;

        for (int row = 0; row < getRows(); row++) {
            if (row == halfSize) {
                leftCopyMat = c21;
                rightCopyMat = c22;
                rowShift = halfSize;
            }
            this.content[row] = Arrays.copyOf(leftCopyMat.content[row
                    - rowShift], getCols());

            for (int col = halfSize; col < getCols(); col++) {
                this.content[row][col] = rightCopyMat.content[row - rowShift][col
                        - halfSize];
            }
        }
    }

    @Override
    public void strassenMultThisWithInto(Matrix matrix, Matrix result) {
        strassenMultThisWithInto(matrix, result, WRITE_BY_ROW);
    }

    @Override
    public void prlStrassenMultThisWithInto(Matrix matrix, Matrix result) {
        prlStrassenMultThisWithInto(matrix, result, WRITE_BY_ROW);
    }

    @Override
    public void winogradMultThisWithInto(Matrix matrix, Matrix result) {
        winogradMultThisWithInto(matrix, result, WRITE_BY_ROW);
    }

    @Override
    public double getMaxNorm() {
        double result = 0;

        for (int row = 0; row < getRows(); ++row) {
            for (int col = 0; col < getCols(); ++col) {
                double entry = Math.abs(content[row][col]);
                if (entry > result) {
                    result = entry;
                }
            }
        }

        return result;
    }

    @Override
    public double get2Norm() {
        double result = 0;

        for (int row = 0; row < getRows(); ++row) {
            for (int col = 0; col < getCols(); ++col) {
                double entry = content[row][col];
                result += entry * entry;
            }
        }

        return Math.sqrt(result);
    }
}
