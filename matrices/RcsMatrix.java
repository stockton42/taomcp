package matrices;

public class RcsMatrix extends Matrix {

    private CcsMatrix content;

    private RcsMatrix(Matrix matrix) {
        this.content = getCcsMatrix(matrix);
    }

    public RcsMatrix(int rows, int cols, int size) {
        this.content = new CcsMatrix(cols, rows, size);
    }

    private CcsMatrix getCcsMatrix(Matrix matrix) {
        if (matrix instanceof CcsMatrix) {
            return (CcsMatrix) matrix;
        } else {
            throw new IllegalArgumentException();
        }
    }

    @Override
    public Matrix clone() {
        return new RcsMatrix(content.clone());
    }

    @Override
    protected Matrix multWith(Matrix matrix) {
        return new RcsMatrix(getCcsMatrix(matrix.clone()).multWith(content));
    }

    @Override
    protected Matrix prlMultWith(Matrix matrix) {
        return new RcsMatrix(getCcsMatrix(matrix.clone()).prlMultWith(content));
    }

    @Override
    public double get(int row, int col) {
        return content.get(col, row);
    }

    @Override
    public void put(double val, int row, int col) {
        content.put(val, col, row);
    }

    @Override
    public void del(int row, int col) {
        content.del(col, row);
    }

    @Override
    public int getRows() {
        return content.getCols();
    }

    @Override
    public int getCols() {
        return content.getRows();
    }

    @Override
    public Matrix getPart(int row1, int col1, int row2, int col2) {
        RcsMatrix result = new RcsMatrix(
                content.getPart(col1, row1, col2, row2));

        return result;
    }

    @Override
    public Matrix getNewInstance(int rows, int cols) {
        return new RcsMatrix(content.getNewInstance(cols, rows));
    }

    @Override
    public boolean isNonNegative() {
        return content.isNonNegative();
    }

    @Override
    protected void setNegativeEntriesToZero(boolean showModifications) {
        content.setNegativeEntriesToZero(showModifications);
    }

    @Override
    public void add(Matrix mat) {
        content.add(mat);
    }

    @Override
    public void stabilizeRowsTo(double stabilizeRowsTo) {
        content.stabilizeColsTo(stabilizeRowsTo);
    }

    @Override
    public void sub(Matrix mat) {
        content.sub(mat);
    }

    @Override
    protected Matrix strassenMultThisWith(Matrix matrix) {
        return new RcsMatrix(getCcsMatrix(matrix.clone()).strassenMultThisWith(
                content));
    }

    @Override
    protected void pool(Matrix upLeft, Matrix upRight, Matrix downLeft,
            Matrix downRight) {
        CcsMatrix upLeftCcs = getCcsMatrix(upLeft);
        CcsMatrix upRightCcs = getCcsMatrix(upRight);
        CcsMatrix downLeftCcs = getCcsMatrix(downLeft);
        CcsMatrix downRightCcs = getCcsMatrix(downRight);

        content.pool(upLeftCcs, upRightCcs, downLeftCcs, downRightCcs);
    }

    @Override
    protected Matrix prlStrassenMultThisWith(Matrix matrix) {
        return new RcsMatrix(getCcsMatrix(matrix.clone())
                .prlStrassenMultThisWith(content));
    }

    @Override
    protected Matrix winogradMultThisWith(Matrix matrix) {
        return new RcsMatrix(getCcsMatrix(matrix.clone()).winogradMultThisWith(
                content));
    }

    /**
     * Only works if norm(A) == norm(transpose(A)).
     */
    @Override
    public double getNorm(MatrixNorm norm) {
        return content.getNorm(norm);
    }

    @Override
    public void stabilizeColsTo(double stabilizeColsTo) {
        content.stabilizeRowsTo(stabilizeColsTo);
    }

}
