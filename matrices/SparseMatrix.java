package matrices;

import java.util.Arrays;

import tools.ArrayHelper;

/**
 * Is a model for a matrix containing only a few double values. Uses the
 * Compressed Column Storage (CCS) algorithm. Best performance is provided if
 * less than 0,1 % of all entries are used. Do NOT use access methods with
 * multiple threads. This matrix is not thread-safe!
 * 
 * @author Michael Stock
 */
public class SparseMatrix extends Matrix {

    private int rows;
    private int cols;

    private double[] val;
    private int[] row_idx;
    private int[] col_ptr;
    private int nextValIndex;
    private int size;

    private static final int ARRAY_MULT_FACTOR = 2;
    private static final int NO_POSITION = -1;
    private static final int DEFAULT_VALUE = 0;
    private static final int NUMBER_OF_THREADS = 8;

    // write complete columns first in multiplication
    private static final boolean WRITE_BY_ROW = false;

    // increases speed when the matrix is not too sparse
    private static final boolean BINARY_ACCESS = true;

    public SparseMatrix(int rows, int cols, int initNumberOfVals) {
        initWith(rows, cols, initNumberOfVals);
    }

    public SparseMatrix(int rows, int cols) {
        initWith(rows, cols, rows + cols);
    }

    private void initWith(int rows, int cols, int initNumberOfVals) {
        if (initNumberOfVals <= 0) {
            throw new IllegalArgumentException();
        }
        this.rows = rows;
        this.cols = cols;
        nextValIndex = 0;
        size = initNumberOfVals;
        initArraysWithSize(size, cols);
    }

    private void initArraysWithSize(int size, int cols) {
        val = new double[size];
        row_idx = new int[size];
        col_ptr = new int[cols + 1];
    }

    public SparseMatrix(Matrix mat) {
        this.rows = mat.getRows();
        this.cols = mat.getCols();
        nextValIndex = 0;
        size = rows + cols;
        initArraysWithSize(size, cols);
        int entryCount = 0;
        for (int c = 0; c < cols; ++c) {
            this.col_ptr[c] = entryCount;
            for (int r = 0; r < rows; ++r) {
                double entry = mat.get(r, c);
                if (entry != DEFAULT_VALUE) {
                    setLastEntryAt(entry, r);
                    entryCount++;
                }
            }
        }
        this.col_ptr[cols] = entryCount;
    }

    private SparseMatrix(SparseMatrix[] matrices) {
        this.rows = matrices[0].rows;
        this.cols = matrices[0].cols;
        size = 0;
        for (int matIndex = 0; matIndex < matrices.length; ++matIndex) {
            size += matrices[matIndex].nextValIndex;
        }
        initArraysWithSize(size, cols);
        nextValIndex = 0;
        int entryCount = 0;

        for (int col = 0; col < cols; ++col) {
            this.col_ptr[col] = entryCount;
            int matIndex = col % matrices.length;
            for (int row = 0; row < rows; ++row) {
                double entry = matrices[matIndex].get(row, col);
                if (entry != DEFAULT_VALUE) {
                    setLastEntryAt(entry, row);
                    entryCount++;
                }
            }
        }
        this.col_ptr[cols] = entryCount;
    }

    public static void main(String[] args) {
        accessTest();
        delTest();
    }

    public static void accessTest() {
        int rows = 1_000;
        int cols = 1_000;
        double fill_percentage = 0.001;

        int vals = (int) ((rows * cols / 100) * fill_percentage);
        long time = System.currentTimeMillis();
        ArrayMatrix am = new ArrayMatrix(rows, cols);
        // SparseMatrix sm = new SparseMatrix(rows, cols, vals);
        for (int j = 0; j < vals; ++j) {
            am.put(1 + ((int) (Math.random() * vals)),
                    (int) (Math.random() * rows), (int) (Math.random() * cols));
        }

        SparseMatrix sm = new SparseMatrix(am);
        if (!sm.equals(am)) {
            throw new IllegalStateException("MATRICES ARE NOT EQUAL!");
        }
        time = System.currentTimeMillis() - time;
        System.out.println("SET UP: " + time + " ms");

        int runs = 1;
        time = System.currentTimeMillis();
        for (int run = 0; run < runs; ++run) {
            for (int c = 0; c < cols; ++c) {
                for (int r = 0; r < rows; ++r) {
                    sm.getPosition(r, c);
                }
            }
        }
        time = System.currentTimeMillis() - time;
        System.out.println("ELAPSED: " + time + " ms");
    }

    public static void delTest() {
        SparseMatrix test = new SparseMatrix(3, 3, 2);
        test.printStatus();
        test.put(1, 0, 0);
        test.printStatus();
        test.put(2, 1, 2);
        test.printStatus();
        test.put(3, 2, 1);
        test.printStatus();
        System.out.println(test);
        test.del(2, 1);
        test.printStatus();
        test.del(1, 2);
        test.printStatus();
        test.del(0, 0);
        test.printStatus();
        System.out.println(test);
    }

    private int getPosition(int row, int col) {
        int startPos = col_ptr[col];
        int maxIndex = col_ptr[col + 1] - 1;

        if (startPos > maxIndex) {
            return NO_POSITION;
        } else if (startPos == maxIndex) {
            if (row_idx[startPos] == row) {
                return startPos;
            } else {
                return NO_POSITION;
            }
        }

        if (!BINARY_ACCESS) {
            for (int i = startPos; i <= maxIndex; ++i) {
                if (row_idx[i] == row) {
                    return i;
                }
            }
        } else {
            // binary search: row_idx is in the right order
            while (startPos <= maxIndex) {
                int i = (maxIndex - startPos) / 2 + startPos;
                if (row_idx[i] > row) {
                    maxIndex = i - 1;
                } else if (row_idx[i] < row) {
                    startPos = i + 1;
                } else {
                    return i;
                }
            }
        }

        return NO_POSITION;
    }

    @Override
    public double get(int row, int col) {
        if (!isValidEntryLocation(row, col))
            return DEFAULT_VALUE;

        int position = getPosition(row, col);
        if (position != NO_POSITION) {
            return val[position];
        } else {
            return DEFAULT_VALUE;
        }
    }

    private void setLastEntryAt(double val, int row) {
        // no entry position test

        if (nextValIndex == size) {
            size *= ARRAY_MULT_FACTOR;
            this.val = Arrays.copyOf(this.val, size);
            this.row_idx = Arrays.copyOf(this.row_idx, size);
        }

        this.val[nextValIndex] = val;
        this.row_idx[nextValIndex] = row;

        nextValIndex++;
    }

    @Override
    public void put(double val, int row, int col) {
        if (!isValidEntryLocation(row, col))
            return;

        int position = getPosition(row, col);
        if (position == NO_POSITION && val != DEFAULT_VALUE) {
            // get next larger row or column index
            int minNewIndex = col_ptr[col];
            int maxNewIndex = col_ptr[col + 1] - 1;

            // enlarge storage if necessary
            if (nextValIndex == size) {
                size *= ARRAY_MULT_FACTOR;
                this.val = Arrays.copyOf(this.val, size);
                this.row_idx = Arrays.copyOf(this.row_idx, size);
            }

            int insertPosition = -1;
            if (maxNewIndex == -1) {
                // empty matrix
                insertPosition = 0;
            } else if (minNewIndex > maxNewIndex) {
                // empty column
                insertPosition = minNewIndex;
            } else if (minNewIndex == maxNewIndex) {
                // one entry
                if (row_idx[minNewIndex] < row) {
                    insertPosition = minNewIndex + 1;
                } else {
                    insertPosition = minNewIndex;
                }
            } else if (row_idx[minNewIndex] > row) {
                // first entry in column
                insertPosition = minNewIndex;
            } else if (row_idx[maxNewIndex] < row) {
                // last entry in column
                insertPosition = maxNewIndex + 1;
            } else {
                // middle entry in column
                if (!BINARY_ACCESS) {
                    for (int i = minNewIndex; i < maxNewIndex; ++i) {
                        if (row_idx[i] < row && row < row_idx[i + 1]) {
                            insertPosition = i + 1;
                            break;
                        }
                    }
                } else {
                    // binary search: row_idx is in the right order
                    while (insertPosition == -1) {
                        int i = (maxNewIndex - minNewIndex) / 2 + minNewIndex;
                        if (row_idx[i] < row && row < row_idx[i + 1]) {
                            insertPosition = i + 1;
                        } else if (row_idx[i] < row) {
                            minNewIndex = i;
                        } else if (row_idx[i + 1] > row) {
                            maxNewIndex = i;
                        } else {
                            throw new IllegalStateException(
                                    "WRONG ORDER IN row_idx!");
                        }
                    }
                }
            }

            for (int i = col + 1; i < this.getCols() + 1; ++i) {
                col_ptr[i]++;
            }

            this.val = ArrayHelper.shift(this.val, insertPosition,
                    nextValIndex + 1, false, true);
            this.row_idx = ArrayHelper.shift(this.row_idx, insertPosition,
                    nextValIndex + 1, false, true);

            this.val[insertPosition] = val;
            row_idx[insertPosition] = row;

            nextValIndex++;
        } else if (position != NO_POSITION && val == DEFAULT_VALUE) {
            if (val == DEFAULT_VALUE) {
                del(row, col);
            }
        } else if (position != NO_POSITION && val != DEFAULT_VALUE) {
            // overwrite
            this.val[position] = val;
        } else { // position == NO_POSITION && val == DEFAULT_VALUE
            // ignore
        }
    }

    @Override
    public void del(int row, int col) {
        if (!isValidEntryLocation(row, col))
            return;

        int position = getPosition(row, col);
        if (position != NO_POSITION) {
            this.val = ArrayHelper.shift(this.val, position, nextValIndex + 1,
                    false, false);
            this.row_idx = ArrayHelper.shift(this.row_idx, position,
                    nextValIndex + 1, false, false);
            nextValIndex--;
            for (int i = col + 1; i < this.getCols() + 1; ++i) {
                col_ptr[i]--;
            }
        }
    }

    @Override
    public int getRows() {
        return rows;
    }

    @Override
    public int getCols() {
        return cols;
    }

    public void printStatus() {
        System.out.println("---STATUS---");
        System.out.println("LENGTH:\t\t " + nextValIndex);
        System.out.println("RESERVED SPACE:\t " + val.length);
        System.out.println("COL_PTR:\t " + ArrayHelper.toString(col_ptr));
        System.out.println("ROW_IDX:\t " + ArrayHelper.toString(row_idx));
        System.out.println("VALS:\t\t " + ArrayHelper.toString(this.val));
    }

    @Override
    public Matrix multWith(Matrix matrix) {
        if (!multPossible(matrix)) {
            throw new IllegalArgumentException();
        }

        SparseMatrix result = new SparseMatrix(this.getRows(),
                matrix.getCols(), this.getRows() + matrix.getCols());

        multThisWithInto(matrix, result, WRITE_BY_ROW);

        return result;
    }

    @Override
    public Matrix prlMultWith(Matrix matrix) {
        if (!multPossible(matrix)) {
            throw new IllegalArgumentException();
        }

        SparseMatrix[] matrices = new SparseMatrix[NUMBER_OF_THREADS];
        for (int i = 0; i < NUMBER_OF_THREADS; ++i) {
            matrices[i] = new SparseMatrix(this.getRows(), matrix.getCols(),
                    this.getRows() + matrix.getCols());
        }

        prlMultThisWithInto(matrix, matrices, NUMBER_OF_THREADS, WRITE_BY_ROW);

        SparseMatrix result = new SparseMatrix(matrices);

        return result;
    }

    @Override
    public Matrix getNewInstance(int rows, int cols) {
        return new SparseMatrix(rows, cols);
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

        SparseMatrix result = new SparseMatrix(row2 - row1 + 1, col2 - col1 + 1);

        if (row1 < 0 || col1 < 0 || row2 < 0 || col2 < 0
                || col1 >= this.getCols() || row1 >= this.getRows()) {
            // invalid arguments: return empty matrix
            return result;
        }

        row2 = Math.min(this.getRows() - 1, row2);
        col2 = Math.min(this.getCols() - 1, col2);

        int entryCount = 0;
        for (int col = col1; col <= col2; ++col) {
            result.col_ptr[col - col1] = entryCount;
            for (int index = col_ptr[col]; index < col_ptr[col + 1]; ++index) {
                if (row1 <= row_idx[index] && row_idx[index] <= row2
                        && val[index] != DEFAULT_VALUE) {
                    result.setLastEntryAt(val[index], row_idx[index] - row1);
                    entryCount++;
                }
            }
        }
        result.col_ptr[col2 - col1 + 1] = entryCount;

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
            SparseMatrix temp = new SparseMatrix(this.getRows(), this.getCols());

            int entryCount = 0;
            double val;
            int nextEntryIndex;

            for (int col = 0; col < this.getCols(); ++col) {
                temp.col_ptr[col] = entryCount;
                nextEntryIndex = nextEntryIndexForColumn(col,
                        this.col_ptr[col] - 1);
                for (int row = 0; row < this.getRows(); ++row) {
                    if (nextEntryIndex != -1 && row_idx[nextEntryIndex] == row) {
                        val = this.val[nextEntryIndex];
                        nextEntryIndex = nextEntryIndexForColumn(col,
                                nextEntryIndex);
                    } else {
                        val = 0;
                    }

                    if (add) {
                        val += mat.get(row, col);
                        if (val != DEFAULT_VALUE) {
                            temp.setLastEntryAt(val, row);
                            entryCount++;
                        }
                    } else {
                        val -= mat.get(row, col);
                        if (val != DEFAULT_VALUE) {
                            temp.setLastEntryAt(val, row);
                            entryCount++;
                        }
                    }
                }
            }
            temp.col_ptr[temp.getCols()] = entryCount;

            this.col_ptr = temp.col_ptr;
            this.val = temp.val;
            this.row_idx = temp.row_idx;
            this.size = temp.size;
            this.nextValIndex = temp.nextValIndex;
        } else {
            throw new IllegalArgumentException();
        }
    }

    private int nextEntryIndexForColumn(int col, int entryIndex) {
        int idx = entryIndex + 1;
        if (idx < col_ptr[col + 1]) {
            return idx;
        } else {
            return -1;
        }
    }

    @Override
    public Matrix clone() {
        SparseMatrix clone = new SparseMatrix(rows, cols, Math.max(size, 1));

        clone.col_ptr = Arrays.copyOf(this.col_ptr, this.col_ptr.length);
        clone.row_idx = Arrays.copyOf(this.row_idx, this.row_idx.length);
        clone.val = Arrays.copyOf(this.val, this.val.length);

        clone.nextValIndex = this.nextValIndex;

        return clone;
    }

    @Override
    public void pool(Matrix upLeft, Matrix upRight, Matrix downLeft,
            Matrix downRight) {
        if (!poolPossible(upLeft, upRight, downLeft, downRight)) {
            throw new IllegalArgumentException();
        }
        SparseMatrix c11 = (SparseMatrix) upLeft;
        SparseMatrix c12 = (SparseMatrix) upRight;
        SparseMatrix c21 = (SparseMatrix) downLeft;
        SparseMatrix c22 = (SparseMatrix) downRight;

        int halfSize = upLeft.getRows();

        size = c11.size + c12.size + c21.size + c22.size;
        initArraysWithSize(size, getCols());
        nextValIndex = 0;
        int entryCount = 0;

        for (int col = 0; col < Math.min(halfSize, getCols()); ++col) {
            col_ptr[col] = entryCount;
            entryCount = copyColFrom(c11, entryCount, col, 0, 0);
            entryCount = copyColFrom(c21, entryCount, col, halfSize, 0);
        }
        for (int col = halfSize; col < getCols(); ++col) {
            col_ptr[col] = entryCount;
            entryCount = copyColFrom(c12, entryCount, col, 0, halfSize);
            entryCount = copyColFrom(c22, entryCount, col, halfSize, halfSize);
        }

        col_ptr[getCols()] = entryCount;
    }

    private int copyColFrom(SparseMatrix mat, int entryCount, int col,
            int rowShift, int colShift) {
        for (int index = mat.col_ptr[col - colShift]; index < mat.col_ptr[col
                + 1 - colShift]; ++index) {
            int row = mat.row_idx[index] + rowShift;
            if (row < getRows()) {
                setLastEntryAt(mat.val[index], row);
                entryCount++;
            }
        }
        return entryCount;
    }

    @Override
    public Matrix strassenMultThisWith(Matrix matrix) {
        SparseMatrix result = new SparseMatrix(this.getRows(),
                matrix.getCols(), this.getRows() + matrix.getCols());

        strassenMultThisWithInto(matrix, result, WRITE_BY_ROW);

        return result;
    }

    @Override
    public Matrix prlStrassenMultThisWith(Matrix matrix) {
        SparseMatrix result = new SparseMatrix(this.getRows(),
                matrix.getCols(), this.getRows() + matrix.getCols());

        prlStrassenMultThisWithInto(matrix, result, WRITE_BY_ROW);

        return result;
    }

    @Override
    public Matrix winogradMultThisWith(Matrix matrix) {
        SparseMatrix result = new SparseMatrix(this.getRows(),
                matrix.getCols(), this.getRows() + matrix.getCols());

        winogradMultThisWithInto(matrix, result, WRITE_BY_ROW);

        return result;
    }

    @Override
    public double getNorm(MatrixNorm norm) {
        switch (norm) {
        case MAX_NORM:
            return getMaxNorm();
        case TWO_NORM:
            return get2Norm();
        default:
            throw new IllegalArgumentException();
        }
    }

    private double getMaxNorm() {
        double result = 0;

        for (int index = 0; index < nextValIndex; ++index) {
            double entry = Math.abs(val[index]);
            if (entry > result) {
                result = entry;
            }
        }

        return result;
    }

    private double get2Norm() {
        double result = 0;

        for (int index = 0; index < nextValIndex; ++index) {
            double entry = val[index];
            result += entry * entry;
        }

        return Math.sqrt(result);
    }

    @Override
    public void stabilizeRowsTo(double stabilizeRowsTo) {
        double[] rowSums = new double[getRows()];

        for (int col = 0; col < getCols(); ++col) {
            for (int colIndex = col_ptr[col]; colIndex < col_ptr[col + 1]; ++colIndex) {
                rowSums[row_idx[colIndex]] += val[colIndex];
            }
        }

        for (int row = 0; row < getRows(); ++row) {
            if (rowSums[row] == 0) {
                rowSums[row] = 1;
            }
        }

        for (int col = 0; col < getCols(); ++col) {
            for (int colIndex = col_ptr[col]; colIndex < col_ptr[col + 1]; ++colIndex) {
                val[colIndex] *= stabilizeRowsTo / rowSums[row_idx[colIndex]];
            }
        }
    }

    @Override
    public boolean isNotNegative() {
        for (int index = 0; index < nextValIndex; ++index) {
            if (val[index] < 0) {
                return false;
            }
        }

        return true;
    }
}
