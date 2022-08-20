import java.io.*;
import java.nio.file.Files;
import java.util.*;

import static javax.swing.UIManager.put;




public class SmithWaterman {
    // Map of amino acids to corresponding index in the BLOSUM62 scoring matrix
    private static final Map<Character, Integer> AMINO_TO_INDEX = new HashMap<>(){
        {put('A', 0);  put('R', 1);  put('N', 2);  put('D', 3);  put('C', 4);
         put('Q', 5);  put('E', 6);  put('G', 7);  put('H', 8);  put('I', 9);
         put('L', 10); put('K', 11); put('M', 12); put('F', 13); put('P', 14);
         put('S', 15); put('T', 16); put('W', 17); put('Y', 18); put('V', 19);
         put('a', 0);  put('r', 1);  put('n', 2);  put('d', 3);  put('c', 4);
         put('q', 5);  put('e', 6);  put('g', 7);  put('h', 8);  put('i', 9);
         put('l', 10); put('k', 11); put('m', 12); put('f', 13); put('p', 14);
         put('s', 15); put('t', 16); put('w', 17); put('y', 18); put('v', 19);}};

    // BLOSUM62 scoring matrix
    private static final int[][] BLOSUM62 = {
        {4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0},
        {-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3},
        {-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3},
        {-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3},
        {0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1},
        {-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2},
        {-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2},
        {0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3},
        {-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3},
        {-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3},
        {-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1},
        {-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2},
        {-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1},
        {-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1},
        {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2},
        {1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2},
        {0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0},
        {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3},
        {-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1},
        {0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4}} ;

    // Linear gap cost (used when a '-' is inserted into the alignment)
    private static final int GAP_COST = -4;

    // Input values
    private String s1id;
    private String s2id;
    private String s1;
    private String s2;
    private int n;

    // Scoring matrix
    private int[][] score;

    // value and location of max value in scoring matrix
    // the max value corresponds to the start of the optimal alignment
    private int max_val;
    private int max_i;
    private int max_j;


    public SmithWaterman(String s1id, String s2id, String s1, String s2, int n) {
        this.s1id = s1id;
        this.s2id = s2id;
        this.s1 = s1;
        this.s2 = s2;
        this.n = n;
        this.score = new int[s1.length() + 1][s2.length() + 1];
        fillScore();
        optimizeAlignment();
    }

    public void fillScore() {
        // fill first column with 0s
        for (int i = 0; i < score.length; i++) {
            score[i][0] = 0;
        }
        // fill first row with 0s
        for (int j = 0; j < score[0].length; j++) {
            score[0][j] = 0;
        }

        // fill remainder of score matrix
        for (int i = 1; i < score.length; i++) {
            for (int j = 1; j < score[0].length; j++) {
                List<Integer> options = new ArrayList<>();
                options.add(score[i - 1][j - 1] + BLOSUM62[AMINO_TO_INDEX.get(s1.charAt(i - 1))][AMINO_TO_INDEX.get(s2.charAt(j - 1))]);
                options.add(score[i - 1][j] + GAP_COST);
                options.add(score[i][j - 1] + GAP_COST);
                options.add(0);
                score[i][j] = Collections.max(options);
            }
        }
    }

    public void optimizeAlignment() {
        max_val = 0;
        max_i = 0;
        max_j = 0;

        // iterate through the array updating max when a larger value is found
        for (int i = 1; i < s1.length(); i++) {
            for (int j = 1; j < s2.length(); j++) {
                if (score[i][j] > max_val) {
                    max_val = score[i][j];
                    max_i = i;
                    max_j = j;
                }
            }
        }
    }

    // return int representing score
    public int score(){
        return max_val;
    }

    // return double representing p-val (iterations will be determined by constructor)
    public double pval(){
        int better = 0;
        for (int k = 0; k < n; k++) {
            Random random = new Random();
            char[] s2chars = s2.toCharArray();
            for (int i = 0; i < s2chars.length; i++) {
                int j = random.nextInt(s2chars.length);
                char temp = s2chars[i];
                s2chars[i] = s2chars[j];
                s2chars[j] = temp;
            }
            String s2random = new String(s2chars);
            SmithWaterman sm = new SmithWaterman("", "", s1, s2random, 0);
            if (sm.score() >= this.score()) {
                better += 1;
            }
        }
        return (better + 1) * 1.0 / (n+1);
    }

    // prints alignment.
    // If both string are less than 15 characters, scoring matrix is printed.
    // if n > 0, p-value is printed
    public void printAlignment() {
        // backtrack for the max_val to find the optimal solution!
        int i = max_i;
        int j = max_j;
        String s1print = "";
        String compare = "";
        String s2print = "";

        // begin by adding the last amino acid to the sequence.
        s1print = s1.charAt(i) + s1print;
        s2print = s2.charAt(j) + s2print;

        // compare that match is s1 and s2 to add to compare string.
        if (s1.charAt(i) == s2.charAt(j)) {
            compare = s1.charAt(i) + compare;
        } else if (BLOSUM62[AMINO_TO_INDEX.get(s1.charAt(i))][AMINO_TO_INDEX.get(s2.charAt(j))] > 0) {
            compare = "+" + compare;
        } else {
            compare = " " + compare;
        }

        // continue to format strings
        while (score[i][j] > 0) {
            if (score[i][j] - GAP_COST == score[i - 1][j]) {
                i--;
                s1print = s1.charAt(i) + s1print;
                s2print = "-" + s2print;
                compare = " " + compare;
            } else if (score[i][j] - GAP_COST == score[i][j - 1]) {
                j--;
                s1print = "-" + s1print;
                s2print = s2.charAt(j) + s2print;
                compare = " " + compare;
            } else {
                i--;
                j--;
                s1print = s1.charAt(i) + s1print;
                s2print = s2.charAt(j) + s2print;
                if (s1.charAt(i) == s2.charAt(j)) {
                    compare = s1.charAt(i) + compare;
                } else if (BLOSUM62[AMINO_TO_INDEX.get(s1.charAt(i))][AMINO_TO_INDEX.get(s2.charAt(j))] > 0) {
                    compare = "+" + compare;
                } else {
                    compare = " " + compare;
                }
            }
        }

        // print some of the basics
        System.out.println("COMPARISON OF " + s1id + " AND " + s2id + "\n");
        System.out.println("Score:" + max_val +"\n");
        System.out.println("Alignment:");

        // print sequence alignment using s1print, compare, s2print.
        // 60 characters per line
        // at the start of line print identifier and location in input string
        String a = s1id + ":\t" + i + "\t";
        String b = "\t\t\t";
        String c = s2id + ":\t" + j + "\t";
        for (int k = 0; k < s1print.length(); k++) {
            if (k != 0 && k % 60 == 0) {
                System.out.println(a);
                System.out.println(b);
                System.out.println(c);
                System.out.println();
                a = s1id + ":\t" + i + "\t";
                b = "\t\t\t";
                c = s2id + ":\t" + j + "\t";
            }
            a += s1print.charAt(k);
            b += compare.charAt(k);
            c += s2print.charAt(k);
            if (s1print.charAt(k) != '-') {
                i++;
            }
            if (s2print.charAt(k) != '-') {
                j++;
            }
        }
        System.out.println(a);
        System.out.println(b);
        System.out.println(c);

        // if both strings are less than 15 letters long, print score matrix
        if (s1.length() < 15 && s2.length() < 15) {
            System.out.println("\nScore Matrix:");
            for (int l = 0; l < score.length; l++) {
                System.out.println(Arrays.toString(score[l]));
            }
        }

        // if n > 0 print p-value
        if (n > 0) {
            System.out.println("\np-value: " + pval());
        }
    }
}