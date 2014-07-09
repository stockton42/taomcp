package tools;

import java.util.LinkedList;
import java.util.List;

public class MathHelper {

    public static void main(String[] args) {
        System.out.println(twoPowersNeededFor(65536));

        System.out.println(twoPowersNeededFor(1023));

        System.out.println(twoPowersNeededFor(1024));
    }

    public static List<Integer> twoPowersNeededFor(int input) {
        if (input < 0) {
            throw new IllegalArgumentException();
        }

        int base = 1;
        int baseIndex = 0;

        while (base <= input) {
            base *= 2;
            baseIndex++;
        }

        List<Integer> result = new LinkedList<Integer>();

        for (int reducer = base / 2; reducer > 0; reducer /= 2) {
            baseIndex--;
            input -= reducer;
            if (input < 0) {
                input += reducer;
            } else {
                result.add(0, baseIndex);
            }
        }

        return result;
    }
}
