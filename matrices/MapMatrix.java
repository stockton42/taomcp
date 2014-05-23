package matrices;

import java.util.HashMap;
import java.util.Map;

/**
 * Is a model for a matrix containing only double values. It uses a map for row
 * access, and then maps the column number to an entry. This implementation is
 * not thread-safe! However, parallel access is possible if the reading/writing
 * sections are disjoint.
 * 
 * @author Michael Stock
 */
public class MapMatrix extends Matrix {

    // calculate complete rows first
    private static final boolean WRITE_BY_ROW = true;

    private static final int DEFAULT_VALUE = 0;
    private static final int NUMBER_OF_THREADS = 8;

    private Map<Integer, Map<Integer, Double>> content;
    private int rows;
    private int cols;

    public MapMatrix(int rows, int cols) {
        initWith(rows, cols);
    }

    private void initWith(int rows, int cols) {
        this.rows = rows;
        this.cols = cols;
        content = new HashMap<Integer, Map<Integer, Double>>();
        for (int i = 0; i < rows; ++i) {
            content.put(i, new HashMap<Integer, Double>());
        }
    }

    public MapMatrix(Matrix mat) {
        initWith(mat.getRows(), mat.getCols());

        double val;
        for (int i = 0; i < rows; ++i) {
            Map<Integer, Double> colMap = content.get(i);
            for (int j = 0; j < cols; ++j) {
                val = mat.get(i, j);
                if (val != DEFAULT_VALUE) {
                    colMap.put(j, val);
                }
            }
        }
    }

    @Override
    public double get(int row, int col) {
        if (!isValidEntryLocation(row, col))
            return DEFAULT_VALUE;

        Map<Integer, Double> colMap = content.get(row);

        Double entry = colMap.get(col);

        if (entry == null) {
            return DEFAULT_VALUE;
        } else {
            return entry;
        }
    }

    @Override
    public void put(double val, int row, int col) {
        if (!isValidEntryLocation(row, col))
            return;

        if (val == DEFAULT_VALUE) {
            del(row, col);
            return;
        }

        Map<Integer, Double> colMap = content.get(row);

        colMap.put(col, val);
    }

    @Override
    public void del(int row, int col) {
        if (!isValidEntryLocation(row, col))
            return;

        Map<Integer, Double> colMap = content.get(row);

        colMap.remove(col);
    }

    @Override
    public int getRows() {
        return rows;
    }

    @Override
    public int getCols() {
        return cols;
    }

    @Override
    public Matrix multWith(Matrix matrix) {
        if (!multPossible(matrix)) {
            throw new IllegalArgumentException();
        }

        MapMatrix result = new MapMatrix(this.getRows(), matrix.getCols());

        multThisWithInto(matrix, result, WRITE_BY_ROW);

        return result;
    }

    @Override
    public Matrix prlMultWith(Matrix matrix) {
        if (!multPossible(matrix)) {
            throw new IllegalArgumentException();
        }

        MapMatrix result = new MapMatrix(this.getRows(), matrix.getCols());
        MapMatrix[] matrices = new MapMatrix[NUMBER_OF_THREADS];
        for (int i = 0; i < NUMBER_OF_THREADS; ++i) {
            matrices[i] = result;
        }

        prlMultThisWithInto(matrix, matrices, NUMBER_OF_THREADS, WRITE_BY_ROW);

        return result;
    }

    @Override
    public Matrix getNewInstance(int rows, int cols) {
        return new MapMatrix(rows, cols);
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

        MapMatrix result = new MapMatrix(row2 - row1 + 1, col2 - col1 + 1);

        if (row1 < 0 || col1 < 0 || row2 < 0 || col2 < 0
                || col1 >= this.getCols() || row1 >= this.getRows()) {
            // invalid arguments: return empty matrix
            return result;
        }

        row2 = Math.min(this.getRows() - 1, row2);
        col2 = Math.min(this.getCols() - 1, col2);

        for (int row = row1; row <= row2; ++row) {
            Map<Integer, Double> colMap = result.content.get(row - row1);
            for (int col = col1; col <= col2; ++col) {
                double value = this.get(row, col);
                if (value != DEFAULT_VALUE) {
                    colMap.put(col - col1, value);
                }
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
                Map<Integer, Double> colMap = this.content.get(row);
                for (int col = 0; col < this.getCols(); ++col) {
                    if (add) {
                        colMap.put(col, this.get(row, col) + mat.get(row, col));
                    } else {
                        colMap.put(col, this.get(row, col) - mat.get(row, col));
                    }
                }
            }
        } else {
            throw new IllegalArgumentException();
        }
    }

    @Override
    public Matrix clone() {
        MapMatrix clone = new MapMatrix(rows, cols);

        for (int row = 0; row < rows; row++) {
            Map<Integer, Double> colMap = clone.content.get(row);
            Map<Integer, Double> thisColMap = this.content.get(row);
            for (int col = 0; col < cols; col++) {
                colMap.put(col, thisColMap.get(col));
            }
        }

        return clone;
    }

    @Override
    public void pool(Matrix upLeft, Matrix upRight, Matrix downLeft,
            Matrix downRight) {
        if (!poolPossible(upLeft, upRight, downLeft, downRight)) {
            throw new IllegalArgumentException();
        }
        MapMatrix c11 = (MapMatrix) upLeft;
        MapMatrix c12 = (MapMatrix) upRight;
        MapMatrix c21 = (MapMatrix) downLeft;
        MapMatrix c22 = (MapMatrix) downRight;

        int halfSize = upLeft.getRows();
        int rowShift = 0;
        MapMatrix leftCopyMat = c11;
        MapMatrix rightCopyMat = c12;

        for (int row = 0; row < getRows(); row++) {
            Map<Integer, Double> colMap = this.content.get(row);
            if (row == halfSize) {
                leftCopyMat = c21;
                rightCopyMat = c22;
                rowShift = halfSize;
            }
            Map<Integer, Double> leftColMap = leftCopyMat.content.get(row
                    - rowShift);
            Map<Integer, Double> rightColMap = rightCopyMat.content.get(row
                    - rowShift);

            for (int col = 0; col < Math.min(halfSize, getCols()); col++) {
                colMap.put(col, leftColMap.get(col));
            }

            for (int col = halfSize; col < getCols(); col++) {
                colMap.put(col, rightColMap.get(col - halfSize));
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
                Double entry = content.get(row).get(col);
                if (entry != null && Math.abs(entry) > result) {
                    result = Math.abs(entry);
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
                Double entry = content.get(row).get(col);
                if (entry != null) {
                    result += entry * entry;
                }
            }
        }

        return Math.sqrt(result);
    }
}
