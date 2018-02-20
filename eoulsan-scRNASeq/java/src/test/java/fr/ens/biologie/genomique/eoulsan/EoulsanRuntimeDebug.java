package fr.ens.biologie.genomique.eoulsan;

import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

import static fr.ens.biologie.genomique.eoulsan.LocalEoulsanRuntime.newEoulsanRuntime;

/**
 * Created by brelurut on 19/02/18.
 */
public class EoulsanRuntimeDebug {
    public static void  initDebugEoulsanRuntime()
        throws IOException, EoulsanException {

        Logger.getLogger(Globals.APP_NAME).getParent().setLevel(Level.OFF);

        if (!EoulsanRuntime.isRuntime()) {

            newEoulsanRuntime(new Settings(true));
        }

        Logger.getLogger(Globals.APP_NAME).setLevel(Level.OFF);
    }
}
