package matrices;

import java.util.Collections;
import java.util.List;

import tools.MathHelper;

public class MatrixPowerer {

    public static final double NO_STABILIZE = 0;

    public static void main(String[] args) {
        double[][] arr = { { 10, 1 }, { 1, 100 } };

        Matrix test = new ArrayMatrix(arr, false);

        test = new SparseMatrix(test);
        // test = new MapMatrix(test);

        System.out.println(test + "\n");

        for (int i = 0; i < 17; ++i) {
            System.out.println("EXPONENT = " + i);
            Matrix powered = stdPower(test, MatrixMultType.NAIVE, i, 1);
            System.out.println(powered);

            System.out.println("LOG_POWER");
            Matrix fastPwd = logPower(test, MatrixMultType.NAIVE, i, 1);
            System.out.println(fastPwd + "\n");
        }
    }

    public static Matrix stdPower(Matrix mat, MatrixMultType matMultType,
            int exponent, double stabilizeRowsTo) {
        if (exponent == 0) {
            Matrix result = mat.getOne();
            if (stabilizeRowsTo != NO_STABILIZE) {
                result.stabilizeRowsTo(stabilizeRowsTo);
            }
            return result;
        }

        Matrix qn = mat.clone();
        if (stabilizeRowsTo != NO_STABILIZE) {
            qn.stabilizeRowsTo(stabilizeRowsTo);
        }

        for (int i = 1; i < exponent; ++i) {
            qn = qn.multWith(mat, matMultType);
            qn.setNegativeEntriesToZero(); // TODO NEW

            if (stabilizeRowsTo != NO_STABILIZE) {
                qn.stabilizeRowsTo(stabilizeRowsTo);
            }
        }

        return qn;
    }

    public static Matrix logPower(Matrix mat, MatrixMultType matMultType,
            int exponent, double stabilizeRowsTo) {
        if (exponent == 0) {
            Matrix result = mat.getOne();
            if (stabilizeRowsTo != NO_STABILIZE) {
                result.stabilizeRowsTo(stabilizeRowsTo);
            }
            return result;
        }

        List<Integer> twoPowers = MathHelper.twoPowersNeededFor(exponent);
        Matrix qn;
        if (twoPowers.contains(0)) {
            qn = mat.clone();
        } else {
            qn = mat.getOne();
        }
        if (stabilizeRowsTo != NO_STABILIZE) {
            qn.stabilizeRowsTo(stabilizeRowsTo);
        }

        int maxPower = Collections.max(twoPowers);
        Matrix temp = mat.clone();
        if (stabilizeRowsTo != NO_STABILIZE) {
            temp.stabilizeRowsTo(stabilizeRowsTo);
        }

        for (int i = 1; i <= maxPower; ++i) {
            temp = temp.multWith(temp, matMultType);
            temp.setNegativeEntriesToZero(); // TODO NEW

            if (stabilizeRowsTo != NO_STABILIZE) {
                temp.stabilizeRowsTo(stabilizeRowsTo);
            }

            if (twoPowers.contains(i)) {
                qn = qn.multWith(temp, matMultType);
                qn.setNegativeEntriesToZero(); // TODO NEW

                if (stabilizeRowsTo != NO_STABILIZE) {
                    qn.stabilizeRowsTo(stabilizeRowsTo);
                }
            }
        }

        return qn;
    }
}
