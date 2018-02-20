package fr.ens.biologie.genomique.eoulsan.modules.preprocessing;



import fr.ens.biologie.genomique.eoulsan.EoulsanRuntimeDebug;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.*;
import java.util.ArrayList;
import java.util.List;




public class TestMerger{

    @Before
    public void setUp() throws Exception {

        EoulsanRuntimeDebug.initDebugEoulsanRuntime();
    }


    @Rule
    public TemporaryFolder folder = new TemporaryFolder();

    @Test
    public void testMerger() {

        final String filesPath = "src/test/files/";
        final String file1Path = filesPath + "expressionstep_output_expression_cell1.tsv";
        final String file2Path = filesPath + "expressionstep_output_expression_cell2.tsv";
        final String file3Path = filesPath + "expressionstep_output_expression_cell3.tsv";
        final String file4Path = filesPath + "expressionstep_output_expression_cell4.tsv";


        final String refPath = filesPath + "mergingResult_ref.tsv";

        try {

            final File merged = folder.newFile("mergingResult.tsv");
            final File ref = new File(refPath);

            List<File> files = new ArrayList<>();
            files.add(new File(file1Path));
            files.add(new File(file2Path));
            files.add(new File(file3Path));
            files.add(new File(file4Path));

            MatrixCreatorModule.merge(files, merged);


            BufferedReader refReader = new BufferedReader(new FileReader(ref));
            BufferedReader mergedReader = new BufferedReader(new FileReader(merged));

            Assert.assertEquals(refReader.readLine(), mergedReader.readLine());

        } catch (Exception e) {
            e.printStackTrace();
        }

    }
}
