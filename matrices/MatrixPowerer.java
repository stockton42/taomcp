package matrices;

import java.util.Collections;
import java.util.List;

import tools.MathHelper;

public class MatrixPowerer {

    public static final double NO_STABILIZE = 0;

    public static void main(String[] args) {
        double[][] arr = { { 9, 1 }, { 1, 99 } };

        Matrix test = new ArrayMatrix(arr, false);

        test = new CrsMatrix(test);
        // test = new MapMatrix(test);

        System.out.println(test + "\n");

        for (int i = 0; i < 17; ++i) {
            System.out.println("EXPONENT = " + i);
            Matrix powered = stdPower(test, MatrixMultType.NAIVE, i, 1.0, true);
            System.out.println(powered);

            System.out.println("LOG_POWER");
            Matrix fastPwd = logPower(test, MatrixMultType.NAIVE, i, 1.0, true);
            System.out.println(fastPwd + "\n");
        }
    }

    public static Matrix stdPower(Matrix mat, MatrixMultType matMultType,
            int exponent, double stabilizeRowsTo,
            boolean setNegativeEntriesToZero) {
        if (exponent == 0) {
            Matrix result = mat.getOne();
            if (setNegativeEntriesToZero)
                result.setNegativeEntriesToZero();
            if (stabilizeRowsTo != NO_STABILIZE)
                result.stabilizeRowsTo(stabilizeRowsTo);

            return result;
        }

        Matrix argumentMatrix = mat;
        if (setNegativeEntriesToZero)
            argumentMatrix.setNegativeEntriesToZero();
        if (stabilizeRowsTo != NO_STABILIZE)
            argumentMatrix.stabilizeRowsTo(stabilizeRowsTo);

        Matrix qn = argumentMatrix.clone();
        if (setNegativeEntriesToZero)
            qn.setNegativeEntriesToZero();
        if (stabilizeRowsTo != NO_STABILIZE)
            qn.stabilizeRowsTo(stabilizeRowsTo);

        for (int i = 1; i < exponent; ++i) {
            qn = qn.multWith(argumentMatrix, matMultType);
            if (setNegativeEntriesToZero)
                qn.setNegativeEntriesToZero();
            if (stabilizeRowsTo != NO_STABILIZE)
                qn.stabilizeRowsTo(stabilizeRowsTo);
        }

        return qn;
    }

    public static Matrix logPower(Matrix mat, MatrixMultType matMultType,
            int exponent, double stabilizeRowsTo,
            boolean setNegativeEntriesToZero) {
        if (exponent == 0) {
            Matrix result = mat.getOne();
            if (setNegativeEntriesToZero) {
                result.setNegativeEntriesToZero();
            }
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
        if (setNegativeEntriesToZero)
            qn.setNegativeEntriesToZero();
        if (stabilizeRowsTo != NO_STABILIZE)
            qn.stabilizeRowsTo(stabilizeRowsTo);

        int maxPower = Collections.max(twoPowers);
        Matrix temp = mat.clone();
        if (setNegativeEntriesToZero)
            temp.setNegativeEntriesToZero();
        if (stabilizeRowsTo != NO_STABILIZE)
            temp.stabilizeRowsTo(stabilizeRowsTo);

        for (int i = 1; i <= maxPower; ++i) {
            temp = temp.multWith(temp, matMultType);

            if (setNegativeEntriesToZero)
                temp.setNegativeEntriesToZero();
            if (stabilizeRowsTo != NO_STABILIZE)
                temp.stabilizeRowsTo(stabilizeRowsTo);

            if (twoPowers.contains(i)) {
                qn = qn.multWith(temp, matMultType);
                if (setNegativeEntriesToZero)
                    qn.setNegativeEntriesToZero();
                if (stabilizeRowsTo != NO_STABILIZE)
                    qn.stabilizeRowsTo(stabilizeRowsTo);
            }
        }

        return qn;
    }
}
