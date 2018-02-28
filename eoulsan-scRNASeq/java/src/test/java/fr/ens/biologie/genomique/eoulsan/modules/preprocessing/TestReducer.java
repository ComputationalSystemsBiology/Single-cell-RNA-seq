package fr.ens.biologie.genomique.eoulsan.modules.preprocessing;

import fr.ens.biologie.genomique.eoulsan.EoulsanRuntimeDebug;
import fr.ens.biologie.genomique.eoulsan.design.io.DesignFormatFinderInputStream;
import org.junit.*;
import org.junit.rules.TemporaryFolder;

import java.io.*;


public class TestReducer {

    @Before
    public void setUp() throws Exception {

        EoulsanRuntimeDebug.initDebugEoulsanRuntime();
    }

    @Rule
    public TemporaryFolder folder = new TemporaryFolder();

    @Test
    public void testReducer() {
            File a = new File("src/test/files/design.txt");

        try {
            File b = folder.newFile("reducedDesign.tsv");
            CellDataExtractorModule.reduceDesign(a, b);

            File ref = new File("src/test/files/reducedDesign_ref.tsv");
            BufferedReader refReader = new BufferedReader(new FileReader(ref));
            BufferedReader redReader = new BufferedReader(new FileReader(b));

            Assert.assertEquals(refReader.readLine(), redReader.readLine());

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    @Test
    public void testDesignFormat() {
        File a = new File("src/test/files/design.txt");

        try {
            // Test design format version
            DesignFormatFinderInputStream formatFinder = new DesignFormatFinderInputStream(new FileInputStream(a));
            int version = formatFinder.getDesignFormatVersion();
            Assert.assertEquals(version, 2);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

}