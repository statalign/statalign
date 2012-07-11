package statalign.postprocess.plugins.contree;

import java.util.BitSet;

public class Utils {

    public static String bitsetToBinary(BitSet set, int n) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n; i++) sb.append(set.get(i) ? "1" : "0");
        return sb.toString();
    }

}
