import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import matrices.ArrayMatrix;
import matrices.MapMatrix;
import matrices.Matrix;
import matrices.SparseMatrix;

public class MultTest {

    public static void main(String[] args) {
        Map<Integer, String> multTypeIds = new HashMap<Integer, String>();
        Map<Integer, String> matrixStorageTypeIds = new HashMap<Integer, String>();
        setUpStringMaps(multTypeIds, matrixStorageTypeIds);

        Set<Integer> multTypesToCalculate = new HashSet<Integer>();
        Set<Integer> matrixStorageTypesToCalculate = new HashSet<Integer>();
        setUpExperiment(multTypesToCalculate, matrixStorageTypesToCalculate);

        Map<Integer, Matrix> results = new HashMap<Integer, Matrix>();

        boolean printResultMatrix = false;
        boolean printDifferences = true;

        int runs = 1;
        int randomNumbers1 = 10000;
        int randomNumbers2 = 10000;
        long time;

        int rows1 = 1000;
        int cols1 = 1000;

        int rows2 = cols1;
        int cols2 = 1000;

        printExperimentParameters(runs, randomNumbers1, randomNumbers2, rows1,
                cols1, rows2, cols2);

        MapMatrix mapM1 = new MapMatrix(rows1, cols1);
        ArrayMatrix arrM1 = new ArrayMatrix(rows1, cols1);
        SparseMatrix sprM1 = new SparseMatrix(rows1, cols1, 10);
        fillMatrices(randomNumbers1, rows1, cols1, mapM1, arrM1, sprM1);

        MapMatrix mapM2 = new MapMatrix(rows2, cols2);
        ArrayMatrix arrM2 = new ArrayMatrix(new double[rows2][cols2], false);
        SparseMatrix sprM2 = new SparseMatrix(rows2, cols2, 10);
        fillMatrices(randomNumbers2, rows2, cols2, mapM2, arrM2, sprM2);

        System.out.println("---\nCALCULATION TIME COMPARISON\n---\n");

        for (int multType : multTypesToCalculate) {

            // double[][] arr1 = { { 0, 1, 2, 3 }, { 1, 2, 3, 4 }, { 2, 3, 4, 5
            // },
            // { 3, 4, 5, 6 } };
            // double[][] arr2 = { { 0, 3, 6, 9 }, { 2, 5, 8, 11 }, { 4, 7, 10,
            // 13
            // },
            // { 6, 9, 12, 15 } };

            System.out.println(multTypeIds.get(multType)
                    + " MATRIX MULTIPLICATION\n---");

            Matrix result = null, mat1, mat2;

            for (int matrixStorageType : matrixStorageTypesToCalculate) {
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

                time = System.currentTimeMillis();
                for (int run = 0; run < runs; ++run) {
                    result = mat1.getNewInstance(rows1, cols2);
                    if (multType == 0) {
                        result = mat1.multWith(mat2);
                    } else if (multType == 1) {
                        result = mat1.prlMultWith(mat2);
                    } else if (multType == 2) {
                        mat1.winogradMultThisWithInto(mat2, result);
                    } else if (multType == 3) {
                        mat1.strassenMultThisWithInto(mat2, result);
                    }
                }
                time = System.currentTimeMillis() - time;
                System.out.println(matrixStorageTypeIds.get(matrixStorageType)
                        + "_MATRIX_TIME:\t " + time + " ms");
                if (printResultMatrix) {
                    System.out.println(result);
                }

                results.put(multType, result);
            }
            System.out.println();
        }

        if (printDifferences) {
            printDifferences(multTypeIds, results);
        }
    }

    private static void setUpExperiment(Set<Integer> multTypesToCalculate,
            Set<Integer> matrixStorageTypesToCalculate) {
        multTypesToCalculate.add(0);
        multTypesToCalculate.add(1);
        multTypesToCalculate.add(2);
        multTypesToCalculate.add(3);
        matrixStorageTypesToCalculate.add(1);
        matrixStorageTypesToCalculate.add(2);
    }

    private static void setUpStringMaps(Map<Integer, String> multTypeIds,
            Map<Integer, String> matrixStorageTypeIds) {
        multTypeIds.put(0, "NAIVE");
        multTypeIds.put(1, "PARALLEL-NAIVE");
        multTypeIds.put(2, "WINOGRAD");
        multTypeIds.put(3, "STRASSEN-NAIVE HYBRID");
        matrixStorageTypeIds.put(0, "MAP");
        matrixStorageTypeIds.put(1, "ARRAY");
        matrixStorageTypeIds.put(2, "SPARSE");
    }

    private static void printExperimentParameters(int runs, int randomNumbers1,
            int randomNumbers2, int rows1, int cols1, int rows2, int cols2) {
        System.out.println("LEFT MATRIX SIZE:\t " + rows1 + " * " + cols1);
        System.out.println("LEFT RANDOM ENTRIES:\t " + randomNumbers1 + "\n");
        System.out.println("RIGHT MATRIX SIZE:\t " + rows2 + " * " + cols2);
        System.out.println("RIGHT RANDOM ENTRIES:\t " + randomNumbers2 + "\n");
        System.out.println("CALCULATION REPETITIONS: " + runs + "\n");
    }

    private static void printDifferences(Map<Integer, String> multTypeIds,
            Map<Integer, Matrix> results) {
        // difference calculation
        System.out.println("---\nDIFFERENCE CALCULATION\n---\n");
        for (Integer multType : results.keySet()) {
            for (Integer multType2 : results.keySet()) {
                if (multType2 > multType) {
                    System.out.println(multTypeIds.get(multType) + " VS "
                            + multTypeIds.get(multType2) + "\n---");
                    Matrix difference = results.get(multType).cloneSub(
                            results.get(multType2));
                    System.out.println("2-NORM_DIFFERENCE:\t "
                            + difference.get2Norm());
                    System.out.println("MAX-NORM_DIFFERENCE:\t "
                            + difference.getMaxNorm());
                    System.out.println();
                }
            }
        }
    }

    public static void fillMatrices(int randomNumbers2, int rows2, int cols2,
            MapMatrix mapM2, ArrayMatrix arrM2, SparseMatrix sprM2) {
        Matrix mat;
        for (int i = 0; i < 3; ++i) {
            if (i == 0) {
                mat = mapM2;
                for (int j = 0; j < randomNumbers2; ++j) {
                    mat.put((int) (Math.random() * randomNumbers2),
                            (int) (Math.random() * rows2),
                            (int) (Math.random() * cols2));
                }
            } else if (i == 1) {
                mat = arrM2;
            } else {
                mat = sprM2;
            }

            if (i > 0) {
                for (int c = 0; c < mapM2.getCols(); ++c) {
                    for (int r = 0; r < mapM2.getRows(); ++r) {
                        mat.put(mapM2.get(r, c), r, c);
                    }
                }
            }

            // System.out.println(mat);
        }
    }
}
