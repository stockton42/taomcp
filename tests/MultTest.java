package tests;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import matrices.ArrayMatrix;
import matrices.CrsMatrix;
import matrices.MapMatrix;
import matrices.Matrix;
import matrices.MatrixMultType;
import matrices.MatrixNorm;
import matrices.MatrixPowerer;

public class MultTest {

    private static final double NO_STABILIZE = MatrixPowerer.NO_STABILIZE;

    public static final double MACHINE_EPSILON = calculateMachineEpsilon();

    private static double calculateMachineEpsilon() {
        int exponent = 0;

        while (1 + Math.pow(2, -exponent) != 1) {
            exponent++;
        }

        return Math.pow(2, -exponent);
    }

    private static final double[][][] matrices = {
            { { 0, 0, 0.5, 0.5 }, { 1.0, 0, 0, 0 },
                    { 1.0 / 3, 1.0 / 3, 1.0 / 3, 0 }, { 0, 0, 1.0, 0 } },
            { { 0.5, 0.5, 0 }, { 0, 0, 1 }, { 0.5, 0.5, 0 } },
            { { 0, 1, 0 }, { 0, 0.5, 0.5 }, { 0.5, 0, 0.5 } },
            { { 0, 0.5, 0.5 }, { 0.5, 0, 0.5 }, { 1, 0, 0 } } };

    public static void main(String[] args) {
        System.out.println("MACHINE EPSILON:\t " + MACHINE_EPSILON);

        Map<Integer, String> matrixStorageTypeIds = setUpStringMaps();

        List<MatrixMultType> multTypesToCalculate = new LinkedList<MatrixMultType>();
        Set<Integer> matrixStorageTypesToCalculate = new HashSet<Integer>();
        setUpExperiment(multTypesToCalculate, matrixStorageTypesToCalculate);

        Map<MatrixMultType, Matrix> results = new HashMap<MatrixMultType, Matrix>();

        // mode == true: calculate left ^ exponent
        // mode == false: calculate left * right
        boolean mode = true;
        boolean useLogPower = true;

        boolean printResultMatrix = true;
        boolean printDifferences = false;
        boolean checkConvergence = false;
        int checkConvergenceMatrixType = 1; // ARRAY, see matrixStorageTypes

        // use non-positive numbers for access to matrices[][][]
        int randomNumbers1 = 0;
        int randomNumbers2 = 0;

        int runs = 1;
        int exponent = 64;

        // set up the dimensions of the matrices
        int rows1 = 4;
        int cols1 = 4;
        int rows2 = cols1;
        int cols2 = 4;
        int initialCrsMatrixSize = 10;

        // set up conditions for stochastic matrices
        double rowSum = 0.0;
        boolean setNegativeEntriesToZero = true;

        // set random number parameters
        long seed1 = 42424242;
        long seed2 = 1354235;

        printExperimentParameters(runs, randomNumbers1, randomNumbers2, rows1,
                cols1, rows2, cols2, rowSum, seed1, seed2);

        // set up left matrix
        MapMatrix mapM1 = new MapMatrix(rows1, cols1);
        ArrayMatrix arrM1 = new ArrayMatrix(rows1, cols1);
        CrsMatrix sprM1 = new CrsMatrix(rows1, cols1, initialCrsMatrixSize);
        fillMatrices(randomNumbers1, mapM1, arrM1, sprM1, rowSum, seed1);

        // set up right matrix
        MapMatrix mapM2 = new MapMatrix(rows2, cols2);
        ArrayMatrix arrM2 = new ArrayMatrix(rows2, cols2);
        CrsMatrix sprM2 = new CrsMatrix(rows2, cols2, 10);
        fillMatrices(randomNumbers2, mapM2, arrM2, sprM2, rowSum, seed2);

        printExperimentInformation(mode, exponent, useLogPower,
                setNegativeEntriesToZero, arrM1, arrM2);

        long time;
        for (MatrixMultType multType : multTypesToCalculate) {
            System.out.println(multType + " MATRIX MULTIPLICATION\n---");
            Matrix result = null, mat1, mat2;

            for (int matrixStorageType : matrixStorageTypesToCalculate) {
                // set matrices to multiply
                if (matrixStorageType == 0) {
                    mat1 = mapM1;
                    mat2 = mapM2;
                } else if (matrixStorageType == 1) {
                    mat1 = arrM1;
                    mat2 = arrM2;
                } else {
                    mat1 = sprM1;
                    mat2 = sprM2;
                }

                // perform calculation
                time = System.currentTimeMillis();
                for (int run = 0; run < runs; ++run) {
                    if (mode) {
                        if (useLogPower) {
                            result = MatrixPowerer.logPower(mat1, multType,
                                    exponent, rowSum, setNegativeEntriesToZero);

                            // for (int i = 0; i < exponent; ++i) {
                            // System.out.println(i + ":");
                            // System.out.println(mat1 + "\n");
                            // mat1 = mat1.multWith(mat1.clone(), multType);
                            // if (rowSum != NO_STABILIZE)
                            // mat1.stabilizeRowsTo(rowSum);
                            // }
                            // System.out.println(exponent + ":");
                            // System.out.println(mat1 + "\n");
                            // result = mat1;
                        } else {
                            result = MatrixPowerer.stdPower(mat1, multType,
                                    exponent, rowSum, setNegativeEntriesToZero);
                        }
                    } else {
                        result = mat1.multWith(mat2, multType);
                        if (setNegativeEntriesToZero)
                            result.setNegativeEntriesToZero();
                        if (rowSum != NO_STABILIZE)
                            result.stabilizeRowsTo(rowSum);
                    }
                }
                time = System.currentTimeMillis() - time;

                // print information, product of non-negative matrices should be
                // non-negative
                System.out.println(matrixStorageTypeIds.get(matrixStorageType)
                        + "_MATRIX_TIME:\t " + time
                        + " ms\t IS NON-NEGATIVE:\t " + result.isNonNegative());
                System.out.println("\t\t\t\t IS POSITIVE:\t\t "
                        + result.isPositive());
                System.out.println("\t\tMINIMAL POSITIVE ENTRY:\t"
                        + result.getMinimalPositiveEntry());
                if (printResultMatrix)
                    System.out.println(result);
                if (setNegativeEntriesToZero)
                    System.out.println();

                // check if all storage types get the same result
                Matrix oldResult = results.get(multType);
                if (oldResult != null && !oldResult.equals(result)) {
                    // throw new IllegalStateException(
                    // "DIFFERENT RESULTS FOR SAME CALCULATION!\nTYPE:\t "
                    // + multType);
                }
                results.put(multType, result);
            }
            System.out.println();
        }

        if (printDifferences)
            printDifferences(results, printResultMatrix);

        if (checkConvergence)
            checkConvergence(results, checkConvergenceMatrixType, mapM1, arrM1,
                    sprM1);

        System.out.println("\n---\nMACHINE EPSILON\n---\n\n" + MACHINE_EPSILON);
    }

    private static void checkConvergence(Map<MatrixMultType, Matrix> results,
            int checkConvergenceMatrixType, MapMatrix mapM1, ArrayMatrix arrM1,
            CrsMatrix sprM1) {
        Matrix mat;
        switch (checkConvergenceMatrixType) {
        case 0:
            mat = mapM1;
            break;
        case 1:
            mat = arrM1;
            break;
        case 2:
            mat = sprM1;
            break;
        default:
            throw new IllegalArgumentException();
        }

        System.out
                .println("\n---\nDIFFERENCE MAX_NORMS TO NEXT MATRIX POWER\n---\n");
        for (MatrixMultType multType : results.keySet()) {
            Matrix matrixPower = results.get(multType);
            Matrix nextMatrixPower = mat.clone();
            nextMatrixPower = nextMatrixPower.multWith(matrixPower, multType);
            nextMatrixPower.sub(matrixPower);
            System.out.print("FOR " + multType + ": \t");
            System.out.println(nextMatrixPower.getNorm(MatrixNorm.MAX_NORM)
                    + "\n");
        }
    }

    private static void printExperimentInformation(boolean mode, int exponent,
            boolean useLogPower, boolean setNegativeEntriesToZero,
            ArrayMatrix arrM1, ArrayMatrix arrM2) {
        System.out.println("LEFT IS NON-NEGATIVE:\t " + arrM1.isNonNegative());
        System.out.println("RIGHT IS NON-NEGATIVE:\t " + arrM2.isNonNegative());
        System.out.println("\nONLY ALLOW NON-NEGATIVE MATRICES:\t "
                + setNegativeEntriesToZero);
        System.out.println("USE FAST MATRIX POWER ALGORITHM:\t " + useLogPower);
        System.out.println("\n---\nCALCULATION TIME COMPARISON");
        if (mode) {
            System.out.println("CALCULATING LEFT ^ " + exponent + " ...");
        } else {
            System.out.println("CALCULATING LEFT * RIGHT MATRIX ...");
        }
        System.out.println("---\n");
    }

    private static void setUpExperiment(
            List<MatrixMultType> multTypesToCalculate,
            Set<Integer> matrixStorageTypesToCalculate) {
        multTypesToCalculate.add(MatrixMultType.NAIVE);
        // multTypesToCalculate.add(MatrixMultType.PARALLEL_NAIVE);
        // multTypesToCalculate.add(MatrixMultType.WINOGRAD);
        // multTypesToCalculate.add(MatrixMultType.STRASSEN_NAIVE_HYBRID);
        // multTypesToCalculate.add(MatrixMultType.PARALLEL_STRASSEN_NAIVE_HYBRID);
        // matrixStorageTypesToCalculate.add(0); // MAP
        matrixStorageTypesToCalculate.add(1); // ARRAY
        // matrixStorageTypesToCalculate.add(2); // CRS
    }

    private static Map<Integer, String> setUpStringMaps() {
        Map<Integer, String> matrixStorageTypeIds = new HashMap<Integer, String>();

        matrixStorageTypeIds.put(0, "MAP");
        matrixStorageTypeIds.put(1, "ARRAY");
        matrixStorageTypeIds.put(2, "CRS");

        return matrixStorageTypeIds;
    }

    private static void printExperimentParameters(int runs, int randomNumbers1,
            int randomNumbers2, int rows1, int cols1, int rows2, int cols2,
            double rowSum, long seed1, long seed2) {
        System.out.println("LEFT MATRIX SIZE:\t " + rows1 + " * " + cols1);
        System.out.println("LEFT RANDOM ENTRIES:\t " + randomNumbers1 + "\n");
        System.out.println("RIGHT MATRIX SIZE:\t " + rows2 + " * " + cols2);
        System.out.println("RIGHT RANDOM ENTRIES:\t " + randomNumbers2 + "\n");
        System.out.println("CALCULATION REPETITIONS: " + runs + "\n");
        System.out.println("ROW SUM:\t\t " + rowSum + "\n");
        System.out.println("FIRST RANDOM NUMBER SEED:\t " + seed1 + "\n");
        System.out.println("SECOND RANDOM NUMBER SEED:\t " + seed2 + "\n");
    }

    private static void printDifferences(Map<MatrixMultType, Matrix> results,
            boolean printResultMatrix) {
        // difference calculation
        System.out.println("---\nDIFFERENCE CALCULATION\n---\n");
        for (MatrixMultType multType : results.keySet()) {
            for (MatrixMultType multType2 : results.keySet()) {
                if (multType2.ordinal() > multType.ordinal()) {
                    System.out.println(multType + " VS " + multType2 + "\n---");
                    Matrix difference = results.get(multType).cloneSub(
                            results.get(multType2));

                    if (printResultMatrix) {
                        System.out.println(difference);
                        System.out.println("---");
                    }

                    for (MatrixNorm norm : MatrixNorm.values()) {
                        System.out.println(norm + " OF DIFFERENCE:\t "
                                + difference.getNorm(norm));
                    }
                    System.out.println();
                }
            }
        }
    }

    public static void fillMatrices(int randomNumbers, MapMatrix mapMat,
            ArrayMatrix arrMat, CrsMatrix sprMat, double rowSum, long seed) {

        Random random = new Random(seed);

        Matrix mat;
        for (int i = 0; i < 3; ++i) {
            if (i == 0) {
                mat = mapMat;
                if (randomNumbers <= 0) {
                    int index = -randomNumbers;
                    for (int row = 0; row < mat.getRows(); ++row) {
                        for (int col = 0; col < mat.getCols(); ++col) {
                            mat.put(matrices[index][row][col], row, col);
                        }
                    }
                } else {
                    for (int j = 0; j < randomNumbers; ++j) {
                        mat.put((int) (random.nextDouble() * randomNumbers),
                                random.nextInt(mat.getRows()),
                                random.nextInt(mat.getCols()));
                    }
                }

                // check sum of all entries in a row
                if (rowSum != NO_STABILIZE) {
                    for (int row = 0; row < mat.getRows(); ++row) {
                        double sum = 0.0;
                        for (int col = 0; col < mat.getCols(); ++col) {
                            sum += mat.get(row, col);
                        }
                        if (sum == 0.0) {
                            mat.put(rowSum, row, row);
                        }
                    }
                    mat.stabilizeRowsTo(rowSum);
                }
            } else if (i == 1) {
                mat = arrMat;
            } else {
                mat = sprMat;
            }

            if (i > 0) {
                copyMatrixFromTo(mapMat, mat);
            }

            // System.out.println(mat);
        }
    }

    private static void copyMatrixFromTo(MapMatrix mapMat, Matrix mat) {
        for (int c = 0; c < mapMat.getCols(); ++c) {
            for (int r = 0; r < mapMat.getRows(); ++r) {
                mat.put(mapMat.get(r, c), r, c);
            }
        }
    }
}
