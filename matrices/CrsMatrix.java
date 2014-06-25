package matrices;

import java.util.Arrays;

import tools.ArrayHelper;

/**
 * Is a model for a matrix containing only a few double values. Uses the
 * Compressed Row Storage (CRS) algorithm. Best performance is provided if less
 * than 0,1 % of all entries are used. Do NOT use access methods with multiple
 * threads. This matrix is not thread-safe!
 * 
 * @author Michael Stock
 */
public class CrsMatrix extends Matrix {

    private int rows;
    private int cols;

    private double[] val;
    private int[] col_idx;
    private int[] row_ptr;
    private int nextValIndex;
    private int size;

    private static final int ARRAY_MULT_FACTOR = 2;
    private static final int NO_POSITION = -1;
    private static final double DEFAULT_VALUE = 0.0;
    private static final int NUMBER_OF_THREADS = 8;

    // write complete rows first in multiplication
    private static final boolean WRITE_BY_ROW = true;

    // increases speed when the matrix is not too sparse
    private static final boolean BINARY_ACCESS = true;

    public CrsMatrix(int rows, int cols, int initNumberOfVals) {
        initWith(rows, cols, initNumberOfVals);
    }

    public CrsMatrix(int rows, int cols) {
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
        initArraysWithSize(size, rows);
    }

    private void initArraysWithSize(int size, int rows) {
        val = new double[size];
        col_idx = new int[size];
        row_ptr = new int[rows + 1];
    }

    public CrsMatrix(Matrix mat) {
        this.rows = mat.getRows();
        this.cols = mat.getCols();
        nextValIndex = 0;
        size = rows + cols;
        initArraysWithSize(size, rows);
        int entryCount = 0;
        for (int r = 0; r < rows; ++r) {
            this.row_ptr[r] = entryCount;
            for (int c = 0; c < cols; ++c) {
                double entry = mat.get(r, c);
                if (entry != DEFAULT_VALUE) {
                    setLastEntryAt(entry, c);
                    entryCount++;
                }
            }
        }
        this.row_ptr[rows] = entryCount;
    }

    private CrsMatrix(CrsMatrix[] matrices) {
        this.rows = matrices[0].rows;
        this.cols = matrices[0].cols;
        size = 0;
        for (int matIndex = 0; matIndex < matrices.length; ++matIndex) {
            size += matrices[matIndex].nextValIndex;
        }
        initArraysWithSize(size, rows);
        nextValIndex = 0;
        int entryCount = 0;

        for (int row = 0; row < rows; ++row) {
            this.row_ptr[row] = entryCount;
            int matIndex = row % matrices.length;
            for (int col = 0; col < cols; ++col) {
                double entry = matrices[matIndex].get(row, col);
                if (entry != DEFAULT_VALUE) {
                    setLastEntryAt(entry, col);
                    entryCount++;
                }
            }
        }
        this.row_ptr[rows] = entryCount;
    }

    public static void main(String[] args) {
        accessTest();
        delTest();
        paperTest();
    }

    private static void paperTest() {
        System.out.println("\n---\nPAPER TEST\n---\n");
        CrsMatrix crs = new CrsMatrix(5, 5, 9);

        crs.put(3.0, 0, 0);
        crs.put(1.0, 0, 2);
        crs.put(2.0, 0, 3);
        crs.put(4.0, 1, 1);
        crs.put(7.0, 2, 1);
        crs.put(5.0, 2, 2);
        crs.put(9.0, 2, 3);
        crs.put(6.0, 4, 3);
        crs.put(5.0, 4, 4);

        System.out.println(crs);
        crs.printStatus();
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

        CrsMatrix sm = new CrsMatrix(am);
        if (!sm.equals(am)) {
            throw new IllegalStateException("MATRICES ARE NOT EQUAL!");
        }
        time = System.currentTimeMillis() - time;
        System.out.println("SET UP: " + time + " ms");

        int runs = 1;
        time = System.currentTimeMillis();
        for (int run = 0; run < runs; ++run) {
            for (int r = 0; r < rows; ++r) {
                for (int c = 0; c < cols; ++c) {
                    sm.getPosition(r, c);
                }
            }
        }
        time = System.currentTimeMillis() - time;
        System.out.println("ELAPSED: " + time + " ms");
    }

    public static void delTest() {
        CrsMatrix test = new CrsMatrix(3, 3, 2);
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
        System.out.println(test);
        test.del(1, 2);
        test.printStatus();
        System.out.println(test);
        test.del(0, 0);
        test.printStatus();
        System.out.println(test);
    }

    private int getPosition(int row, int col) {
        int startPos = row_ptr[row];
        int maxIndex = row_ptr[row + 1] - 1;

        if (startPos > maxIndex) {
            return NO_POSITION;
        } else if (startPos == maxIndex) {
            if (col_idx[startPos] == col) {
                return startPos;
            } else {
                return NO_POSITION;
            }
        }

        if (!BINARY_ACCESS) {
            for (int i = startPos; i <= maxIndex; ++i) {
                if (col_idx[i] == col) {
                    return i;
                }
            }
        } else {
            // binary search: row_idx is in the right order
            while (startPos <= maxIndex) {
                int i = (maxIndex - startPos) / 2 + startPos;
                if (col_idx[i] > col) {
                    maxIndex = i - 1;
                } else if (col_idx[i] < col) {
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

    private void setLastEntryAt(double val, int col) {
        // no entry position test

        if (nextValIndex == size) {
            size *= ARRAY_MULT_FACTOR;
            this.val = Arrays.copyOf(this.val, size);
            this.col_idx = Arrays.copyOf(this.col_idx, size);
        }

        this.val[nextValIndex] = val;
        this.col_idx[nextValIndex] = col;

        nextValIndex++;
    }

    @Override
    public void put(double val, int row, int col) {
        if (!isValidEntryLocation(row, col))
            return;

        int position = getPosition(row, col);
        if (position == NO_POSITION && val != DEFAULT_VALUE) {
            // get next larger row or column index
            int minNewIndex = row_ptr[row];
            int maxNewIndex = row_ptr[row + 1] - 1;

            // enlarge storage if necessary
            if (nextValIndex == size) {
                size *= ARRAY_MULT_FACTOR;
                this.val = Arrays.copyOf(this.val, size);
                this.col_idx = Arrays.copyOf(this.col_idx, size);
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
                if (col_idx[minNewIndex] < col) {
                    insertPosition = minNewIndex + 1;
                } else {
                    insertPosition = minNewIndex;
                }
            } else if (col_idx[minNewIndex] > col) {
                // first entry in column
                insertPosition = minNewIndex;
            } else if (col_idx[maxNewIndex] < col) {
                // last entry in column
                insertPosition = maxNewIndex + 1;
            } else {
                // middle entry in column
                if (!BINARY_ACCESS) {
                    for (int i = minNewIndex; i < maxNewIndex; ++i) {
                        if (col_idx[i] < col && col < col_idx[i + 1]) {
                            insertPosition = i + 1;
                            break;
                        }
                    }
                } else {
                    // binary search: row_idx is in the right order
                    while (insertPosition == -1) {
                        int i = (maxNewIndex - minNewIndex) / 2 + minNewIndex;
                        if (col_idx[i] < col && col < col_idx[i + 1]) {
                            insertPosition = i + 1;
                        } else if (col_idx[i] < col) {
                            minNewIndex = i;
                        } else if (col_idx[i + 1] > col) {
                            maxNewIndex = i;
                        } else {
                            throw new IllegalStateException(
                                    "WRONG ORDER IN col_idx!");
                        }
                    }
                }
            }

            for (int i = row + 1; i < this.getRows() + 1; ++i) {
                row_ptr[i]++;
            }

            this.val = ArrayHelper.shift(this.val, insertPosition,
                    nextValIndex + 1, false, true);
            this.col_idx = ArrayHelper.shift(this.col_idx, insertPosition,
                    nextValIndex + 1, false, true);

            this.val[insertPosition] = val;
            col_idx[insertPosition] = col;

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
            this.val = ArrayHelper.shift(this.val, position, nextValIndex,
                    false, false);
            this.col_idx = ArrayHelper.shift(this.col_idx, position,
                    nextValIndex, false, false);
            nextValIndex--;
            for (int i = row + 1; i < this.getRows() + 1; ++i) {
                row_ptr[i]--;
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
        System.out.println("COL_PTR:\t " + ArrayHelper.toString(row_ptr));
        System.out.println("ROW_IDX:\t " + ArrayHelper.toString(col_idx));
        System.out.println("VALS:\t\t " + ArrayHelper.toString(this.val));
    }

    @Override
    public Matrix multWith(Matrix matrix) {
        if (!multPossible(matrix)) {
            throw new IllegalArgumentException();
        }

        CrsMatrix result = new CrsMatrix(this.getRows(), matrix.getCols(),
                this.getRows() + matrix.getCols());

        multThisWithInto(matrix, result, WRITE_BY_ROW);

        return result;
    }

    @Override
    public Matrix prlMultWith(Matrix matrix) {
        if (!multPossible(matrix)) {
            throw new IllegalArgumentException();
        }

        CrsMatrix[] matrices = new CrsMatrix[NUMBER_OF_THREADS];
        for (int i = 0; i < NUMBER_OF_THREADS; ++i) {
            matrices[i] = new CrsMatrix(this.getRows(), matrix.getCols(),
                    this.getRows() + matrix.getCols());
        }

        prlMultThisWithInto(matrix, matrices, NUMBER_OF_THREADS, WRITE_BY_ROW);

        CrsMatrix result = new CrsMatrix(matrices);

        return result;
    }

    @Override
    public Matrix getNewInstance(int rows, int cols) {
        return new CrsMatrix(rows, cols);
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

        CrsMatrix result = new CrsMatrix(row2 - row1 + 1, col2 - col1 + 1);

        if (row1 < 0 || col1 < 0 || row2 < 0 || col2 < 0
                || col1 >= this.getCols() || row1 >= this.getRows()) {
            // invalid arguments: return empty matrix
            return result;
        }

        row2 = Math.min(this.getRows() - 1, row2);
        col2 = Math.min(this.getCols() - 1, col2);

        int entryCount = 0;
        for (int row = row1; row <= row2; ++row) {
            result.row_ptr[row - row1] = entryCount;
            for (int index = row_ptr[row]; index < row_ptr[row + 1]; ++index) {
                if (col1 <= col_idx[index] && col_idx[index] <= col2
                        && val[index] != DEFAULT_VALUE) {
                    result.setLastEntryAt(val[index], col_idx[index] - col1);
                    entryCount++;
                }
            }
        }
        result.row_ptr[row2 - row1 + 1] = entryCount;

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
            CrsMatrix temp = new CrsMatrix(this.getRows(), this.getCols());

            int entryCount = 0;
            double val;
            int nextEntryIndex;

            for (int row = 0; row < this.getRows(); ++row) {
                temp.row_ptr[row] = entryCount;
                nextEntryIndex = nextEntryIndexForColumn(row,
                        this.row_ptr[row] - 1);
                for (int col = 0; col < this.getCols(); ++col) {
                    if (nextEntryIndex != -1 && col_idx[nextEntryIndex] == col) {
                        val = this.val[nextEntryIndex];
                        nextEntryIndex = nextEntryIndexForColumn(row,
                                nextEntryIndex);
                    } else {
                        val = 0;
                    }

                    if (add) {
                        val += mat.get(row, col);
                        if (val != DEFAULT_VALUE) {
                            temp.setLastEntryAt(val, col);
                            entryCount++;
                        }
                    } else {
                        val -= mat.get(row, col);
                        if (val != DEFAULT_VALUE) {
                            temp.setLastEntryAt(val, col);
                            entryCount++;
                        }
                    }
                }
            }
            temp.row_ptr[temp.getRows()] = entryCount;

            this.row_ptr = temp.row_ptr;
            this.val = temp.val;
            this.col_idx = temp.col_idx;
            this.size = temp.size;
            this.nextValIndex = temp.nextValIndex;
        } else {
            throw new IllegalArgumentException();
        }
    }

    private int nextEntryIndexForColumn(int row, int entryIndex) {
        int idx = entryIndex + 1;
        if (idx < row_ptr[row + 1]) {
            return idx;
        } else {
            return -1;
        }
    }

    @Override
    public Matrix clone() {
        CrsMatrix clone = new CrsMatrix(rows, cols, Math.max(size, 1));

        clone.row_ptr = Arrays.copyOf(this.row_ptr, this.row_ptr.length);
        clone.col_idx = Arrays.copyOf(this.col_idx, this.col_idx.length);
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
        CrsMatrix c11 = (CrsMatrix) upLeft;
        CrsMatrix c12 = (CrsMatrix) upRight;
        CrsMatrix c21 = (CrsMatrix) downLeft;
        CrsMatrix c22 = (CrsMatrix) downRight;

        int halfSize = upLeft.getCols();

        size = c11.size + c12.size + c21.size + c22.size;
        initArraysWithSize(size, getRows());
        nextValIndex = 0;
        int entryCount = 0;

        for (int row = 0; row < Math.min(halfSize, getRows()); ++row) {
            row_ptr[row] = entryCount;
            entryCount = copyRowFrom(c11, entryCount, row, 0, 0);
            entryCount = copyRowFrom(c12, entryCount, row, halfSize, 0);
        }
        for (int row = halfSize; row < getRows(); ++row) {
            row_ptr[row] = entryCount;
            entryCount = copyRowFrom(c21, entryCount, row, 0, halfSize);
            entryCount = copyRowFrom(c22, entryCount, row, halfSize, halfSize);
        }

        row_ptr[getRows()] = entryCount;
    }

    private int copyRowFrom(CrsMatrix mat, int entryCount, int row,
            int colShift, int rowShift) {
        for (int index = mat.row_ptr[row - rowShift]; index < mat.row_ptr[row
                + 1 - rowShift]; ++index) {
            int col = mat.col_idx[index] + colShift;
            if (col < getCols()) {
                setLastEntryAt(mat.val[index], col);
                entryCount++;
            }
        }
        return entryCount;
    }

    @Override
    public Matrix strassenMultThisWith(Matrix matrix) {
        CrsMatrix result = new CrsMatrix(this.getRows(), matrix.getCols(),
                this.getRows() + matrix.getCols());

        strassenMultThisWithInto(matrix, result, WRITE_BY_ROW);

        return result;
    }

    @Override
    public Matrix prlStrassenMultThisWith(Matrix matrix) {
        CrsMatrix result = new CrsMatrix(this.getRows(), matrix.getCols(),
                this.getRows() + matrix.getCols());

        prlStrassenMultThisWithInto(matrix, result, WRITE_BY_ROW);

        return result;
    }

    @Override
    public Matrix winogradMultThisWith(Matrix matrix) {
        CrsMatrix result = new CrsMatrix(this.getRows(), matrix.getCols(),
                this.getRows() + matrix.getCols());

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
        double rowSum;
        for (int row = 0; row < rows; ++row) {
            rowSum = 0;
            for (int colIndex = row_ptr[row]; colIndex < row_ptr[row + 1]; ++colIndex) {
                rowSum += val[colIndex];
            }

            if (rowSum == 0) {
                rowSum = 1;
            }

            for (int colIndex = row_ptr[row]; colIndex < row_ptr[row + 1]; ++colIndex) {
                val[colIndex] *= stabilizeRowsTo / rowSum;
            }
        }
    }

    @Override
    public boolean isNonNegative() {
        for (int index = 0; index < nextValIndex; ++index) {
            if (val[index] < 0) {
                return false;
            }
        }

        return true;
    }

    @Override
    public boolean isPositive() {
        for (int row = 0; row < getRows(); ++row) {
            for (int col = 0; col < getCols(); ++col) {
                if (this.get(row, col) <= 0) {
                    return false;
                }
            }
        }

        return true;
    }

    @Override
    protected void setNegativeEntriesToZero(boolean showModifications) {
        double minValueSetToZero = 0.0;
        for (int index = 0; index < nextValIndex; ++index) {
            if (val[index] < 0) {
                if (showModifications && minValueSetToZero > val[index]) {
                    minValueSetToZero = val[index];
                }
                val[index] = DEFAULT_VALUE;
            }
        }
        if (showModifications && minValueSetToZero < 0) {
            System.out.println("MINIMAL NEGATIVE ENTRY SET TO ZERO: "
                    + minValueSetToZero);
        }
    }

    @Override
    public double getMinimalPositiveEntry() {
        double minimum = 2.0;
        for (int i = 0; i < val.length; ++i) {
            if (minimum > val[i] && val[i] > 0) {
                minimum = val[i];
            }
        }

        return minimum;
    }
}
