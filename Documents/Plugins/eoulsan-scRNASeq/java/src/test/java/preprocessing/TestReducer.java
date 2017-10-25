package preprocessing;

import org.junit.Test;

import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;



public class TestReducer {

    private static void reduceDesign(File file, File output) throws IOException {
        /* Initialise variables */
        String aLine;
        boolean table = false;
        String pattern = "(Columns)";
        Pattern seek = Pattern.compile(pattern);

		/* Read file, copy lines */
        try (BufferedReader in = new BufferedReader(new FileReader(file))) {
            while ((aLine = in.readLine()) != null) {
                if (table) {
                    try (BufferedWriter out = new BufferedWriter(new FileWriter(output, true))) {
                        out.write(aLine);
                        out.newLine();
                    }
                } else {
                    Matcher m = seek.matcher(aLine);
                    if (m.find()) {
                        table = true;
                    }
                }
            }
        }
    }

    @Test
    public void testReducer() {
            File a = new File("../../files/design.txt");
            File b = new File ("../../files/reducedDesign.tsv");
        try {
            reduceDesign(a, b);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }


    }