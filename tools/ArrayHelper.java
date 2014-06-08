package tools;

public class ArrayHelper {

    public static void main(String[] args) {
        double[] toShift = { 1.0, 2.0, 3.0, 4.0, 5.0 };

        toShift = ArrayHelper.shift(toShift, 1, 3, false, true);
        System.out.println(toString(toShift));
    }

    public static double[] shift(double[] toShift, int startPosition,
            int targetPosition, boolean cyclic, boolean rightShift) {
        double left;
        if (cyclic) {
            if (rightShift) {
                left = toShift[targetPosition - 1];
            } else {
                left = toShift[startPosition];
            }
        } else {
            left = 0;
        }
        double swap;
        int index;
        for (int i = 0; i < targetPosition - startPosition; ++i) {
            if (rightShift) {
                index = i + startPosition;
            } else {
                index = targetPosition - 1 - i;
            }
            swap = toShift[index];
            toShift[index] = left;
            left = swap;
        }

        return toShift;
    }

    public static int[] shift(int[] toShift, int startPosition,
            int targetPosition, boolean cyclic, boolean rightShift) {
        int left;
        if (cyclic) {
            if (rightShift) {
                left = toShift[targetPosition - 1];
            } else {
                left = toShift[startPosition];
            }
        } else {
            left = 0;
        }
        int swap;
        int index;
        for (int i = 0; i < targetPosition - startPosition; ++i) {
            if (rightShift) {
                index = i + startPosition;
            } else {
                index = targetPosition - 1 - i;
            }
            swap = toShift[index];
            toShift[index] = left;
            left = swap;
        }

        return toShift;
    }

    public static String toString(int[] array) {
        StringBuilder result = new StringBuilder("{");

        if (array != null && array.length > 0) {
            for (int i = 0; i < array.length - 1; i++) {
                result.append(array[i] + ", ");
            }

            result.append(array[array.length - 1]);
        }

        result.append("}");

        return result.toString();
    }

    public static String toString(double[] array) {
        StringBuilder result = new StringBuilder("{");

        if (array != null && array.length > 0) {
            for (int i = 0; i < array.length - 1; i++) {
                result.append(array[i] + ", ");
            }

            result.append(array[array.length - 1]);
        }

        result.append("}");

        return result.toString();
    }
}
