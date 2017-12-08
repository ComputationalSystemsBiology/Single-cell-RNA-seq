package preprocessing;



import org.junit.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class TestMerger{

    /**
     * Merge two columns files into one matrix
     * @param files list of files to merge
     * @param mergedFile file for writing mergin
     * @throws NullPointerException
     * @throws IOException
     */

    private static void merge(List<File> files, File mergedFile)
            throws NullPointerException, IOException {

			/* Test parameters */
        if (files == null) {
            throw new NullPointerException("files argument cannot be null");
        }

        if (mergedFile == null) {
            throw new NullPointerException("mergedFile argument cannot be null");
        }

			/* Set structures and extract first file */
        File init = files.get(0);

        String headers = init.getName().replace(".tsv", "").split("_")[3];
        List<String> genes = new ArrayList<>();
        Map<String, String> lines = new HashMap<>();

			/* Treat first file initialising structures */
        try (BufferedReader in = new BufferedReader(new FileReader(init))) {
            String aLine;
            String[] parts = new String[2];
            int i = 0;
            while ((aLine = in.readLine()) != null) {
                if (i > 0) {
                    parts = aLine.split("\t");
                    genes.add(parts[0]);
                    lines.put(parts[0], parts[1]);
                }
                i++;
            }
        }

			/* Remove file treated from files list */
        files.remove(0);

			/* Treat other files */
        for (File file : files) {

            headers = headers + "\t" + file.getName().replace(".tsv", "").split("_")[3];
            try (BufferedReader in = new BufferedReader(new FileReader(file))) {
                String aLine;
                String[] parts = new String[2];
                while ((aLine = in.readLine()) != null) {
                    parts = aLine.split("\t");
                    String value = lines.get(parts[0]);
                    value = value + "\t" + parts[1];
                    lines.put(parts[0], value);
                }
            }
        }

			/* Write merged file, excluding non detected genes */
        try (BufferedWriter out = new BufferedWriter(new FileWriter(mergedFile,
                true))) {
            out.write(headers);
            out.newLine();
            for (String key : genes) {
                int sum = 0;
                String[] values = lines.get(key).split("\t");
                for (String value : values) {
                    sum = sum + Integer.parseInt(value);
                }
                if (sum > 0) {
                    out.write(key + "\t" + lines.get(key));
                    out.newLine();
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    @Test
    public void testMerger() {
        File design = new File("../../files/design.txt");

        File merged = new File("../../files/merged.tsv");
        List<File> files = new ArrayList<>();
        files.add(new File ("../../files/expressionstep_output_expression_BRBS99.tsv"));
        files.add(new File ("../../files/expressionstep_output_expression_BRBS100.tsv"));
        files.add(new File ("../../files/expressionstep_output_expression_BRBS114.tsv"));
        files.add(new File ("../../files/expressionstep_output_expression_BRBS122.tsv"));

        try {
            merge(files, merged );
        } catch (Exception e) {
            e.printStackTrace();
        }

    }
}
