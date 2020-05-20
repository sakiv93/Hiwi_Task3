##!/bin/sh
##
## abaquser.sh
## Shell-Skript im Paket ABAQUSER 2.0
##
## Skript zum Aufruf von ABAQUS-Jobs mit UEL- und UMAT-f90-Subroutinen
## Die beiden f90-Dateien werden in einer *.f zusammengefasst.
## ABAQUS-Optionen job= und cpus= werden erkannt, weitere Optionen
## koennen mit opt=irgendeineOption=Wert eingegeben werden.
##
## Ueber ein PYTHON-Skript wird die .inp um die Befehle zum Erzeugen
## der *.fil erweitert.
## Erweiterung der *.odb um die Infos des *.fil, also SDV und diverse
## Knotenwerte ueber ein PYTHON-Skript.
## Ungenutzte ABAQUS-Dateien werden geloescht.
##
## Alle Optionen ueber abaquser help
##
## Skript ausfuehrbarmachen mit chmod +x abaquser.sh
##
## Das Skript ist nicht ausgetestet.
##
## Die Verwendung des Skriptes ist freigegeben
## Aenderungen im Skript sind zu kennzeichnen, im Header zu vermerken
## und dem Autor mitzuteilen.
##
## Stephan Roth
## Stephan.Roth@imfd.tu-freiberg.de
##
## Freiberg, 12.11.2010
##
##########################################################################
##
##########################################################################
## Aenderung 15.12.2010, Stephan Roth:
##   - Paralleles Auslesen des *.fil und Schreiben des *.odb
##     in 2 Prozessen (dazu ABAQUS 6.10.2 notwendig, da dort Python 2.6)
## Aenderung 07.03.2011, Stephan Roth:
##   - Anwendung von ABAQUSER auf GALILEI, Start mit qsub shellskript.sh,
##     in dem abaquser.sh aufgerufen wird
## Aenderung 07.07.2011, Stephan Roth:
##   - Auswahl ob bestehende odb um UEL erweitert werden soll, oder ob sie
##     durch eine odb nur mit UEL ersetzt werden soll
## Aenderung 20.07.2011, Stephan Roth:
##   - vorzeitiges Beenden eines Jobs auf Galilei mit terminate
## Aenderung 18.10.2012, Stephan Roth:
##   - unterbrechen eines Jobs auf Galilei mit suspend
##   - wiederaufnehmen eines Jobs auf Galilei mit resume
## Aenderung 21.07.2011, Stephan Roth:
##   - kopieren aller Fortran-Datein ins Scratch
## Aenderung 02.09.2011, Stephan Roth:
##   - Momentaufnahme eines Jobs auf Galilei mit snapshot
## Aenderung 29.11.2011, Geralf Huetter:
##   - bei Verwendung "opt=oldjob=XYZ" (Restart-Funktion) in Verbindung mit
##     "galilei" werden nun auch die benoetigten Dateien XYZ.* mit in das
##     entsprechende Scratch-Verzeichnis kopiert
## Aenderung 30.11.2011, Geralf Huetter:
##   - bei Verwendung auf galilei werden nun alle *.env-Dateien mitkopiert
## Aenderung 17.10.2011, Stephan Roth:
##   - Vereinigung der Skripte abaquser.sh und abaquser_GalileiBatch.sh
##   - Optionen: env, usr (usr loest umat und uel ab)
##   - Angabe aller Dateien mit Endung und optional Pfadangabe bei Optionen
##     input, usr, env, uelinfo
##   - Angabe mehrerer Dateien in diesen Optionen mit Trennung durch Komma
##   - Auswertung der Aufrufoptionen von fil2odb
## Aenderung 16.01.2012, Stephan Roth:
##   - Multiprocessing von fil2odb: Prozessoranzahl
## Aenderung 20.01.2012, Stephan Roth:
##   - falls Berechnung und fil2odb direkt nacheinander ausgefuehrt werden
##     soll (galilei): kein Zwischenspeichen der (grossen) *.fil und *.odb
##     mehr im home
## Aenderung 20.01.2012, Stephan Roth:
##   - falls Option sharedlib zur Angabe des Pfades fuer ABAQUS shared library
##     libstandardU.so ; Option sharedlib alternativ zu usr (Vorrang sharedlib)
## Aenderung 06.11.2012, Stephan Roth:
##   - modiinp und fil2odb als Userelement-Klasse
## Aenderung 02.04.2013, Stephan Roth:
##   - Vereinfachungen der letzten Aenderung
## Aenderung 03.05.2013, Stephan Roth:
##   - Wahl der Warteschlangen
## Aenderung 04.11.2013, Stephan Roth:
##   - urz-Option fuer URZ HPC-Cluster
## Aenderung 25.04.2014, Stephan Roth:
##   - Loeschen unnuetzer Dateien ueberarbeitet
## Aenderung 15.05.2014, Stephan Roth:
##   - URZ-Cluster Optionen: Anzahl der Chunks (Nodes)
## Aenderung 15.07.2015, Stephan Roth:
##   - Uebergabe des ABAQUS-Aufrufs an Fil2Odb
## Aenderung 04.09.2015, Stephan Roth:
##   - modifiedupdate-Option
## Aenderung 15.02.2018, Stephan Roth:
##   - datacheck
## Aenderung 17.04.2019, Stephan Roth:
##   - mem_free bei Galilei, z.B mem=20G (nicht: mem=20GB !!!)
## Aenderung 23.04.2019, Stephan Roth:
##   - kopieren der *.o ins scratch
## Aenderung 03.05.2019, Stephan Roth:
##   - objectfile-Option
##########################################################################
##
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
number () {
  num=`echo  $1 | awk -v FS="," '{print NF}'`
  echo -e "$num"
}
max () {
  maxvalue=0
  for i in $*
    do
    if [ $i -gt $maxvalue ]
      then
      maxvalue=$i
    fi
  done
  echo -e "$maxvalue"
}
getpathandfile () {
  numa=`echo  $1 | awk -v FS="," '{print NF}'`
  if [ $numa -ge $2 ]
    then
    filestring=`echo  $1 | awk -v FS="," '{print $VAR}'  "VAR=$2"`
    numb=`echo  $filestring | awk -v FS="/" '{print NF}'`
    if [ $numb -eq 1 ]
      then
      path="./"
      file=$filestring
    else
      file=`echo  $filestring | awk -v FS="/" '{print $VAR}' "VAR=$numb"`
      path=`echo  $filestring | awk -v FS="$file" '{print $1}'`
    fi
  else
    echo "$2"" exceeds ""$numa"
    echo "Check call of getpathandfile"
  fi
  if [ "$3" = "p" ]
    then
    echo  "$path"
  elif [ "$3" = "f" ]
    then
    echo  "$file"
  elif [ "$3" = "pf" ]
    then
    echo  "$path""$file"
  fi
}
echo -e "\n"
echo -e "Shellskript zum Start eines ABAQUS-Jobs mit UEL und UMAT"
echo -e "Version:     ABAQUSER 2.0"
echo -e "Infos unter: ABAQUSER.sh help"
date
latestabacall="abq6123"
queuedefaults="default,adde"
if [ "$1" = "help" ]
  then
  echo -e "Eingabeoptionen:"
  echo -e "   ABAQUSER.sh job=jobname           -->" \
          "startet den ABAQUS-Job jobname"
  echo -e "   ABAQUSER.sh input=inputfilename   -->" \
          "Angabe der Inputdatei inputfilename.inp"
  echo -e "   ABAQUSER.sh cpus=anzahl           -->" \
          "entspricht der ABAQUS-Option cpus=anzahl"
  echo -e "   ABAQUSER.sh usr=path1/usersubroutine1.f90,path2/usersubroutine2.f90   -->" \
          "Angabe der Fortran-Dateien mit UserSubRoutinen"
  echo -e "   ABAQUSER.sh sharedlib=path/sharedlibrary.so   -->" \
          "Angabe der ABAQUS shared library libstandardU.so"
  echo -e "   ABAQUSER.sh env=path/abaqus_v6.env-->" \
          "Angabe der ABAQUS Umgebungsdatei"
  echo -e "   ABAQUSER.sh uelinfo=path/uel.info -->" \
          "Angabe der Datei mit UEL Informationen uel.info"
  echo -e "   ABAQUSER.sh opt=abaqusoption=wert -->" \
          "Eingabe weiterer ABAQUS-Optionen, z.B. opt=oldjob=oldjobname"
  echo -e "   ABAQUSER.sh replace               -->" \
          "Loeschen der gleichnamigen job-Daten vor der Berechnung"
  echo -e "   ABAQUSER.sh fil                   -->" \
          "Auslesen des .fil und Erweiterung der .odb"
  echo -e "   ABAQUSER.sh mpfil                 -->" \
          "Auslesen des .fil und Erweiterung der .odb auf 2 Prozessoren"
  echo -e "   ABAQUSER.sh abacall=abqversion    -->" \
          "startet ABAQUS mit Befehl abqversion"
  echo -e "   ABAQUSER.sh noanalysis            -->" \
          "unterdrueckt ABAQUS Rechnung"
  echo -e "   ABAQUSER.sh datacheck               -->" \
          "startet datacheck, keine Rechnung"
  echo -e "   ABAQUSER.sh galilei               -->" \
          "Ausfuehren des Skriptes auf Galilei mit qsub"
  echo -e "   ABAQUSER.sh queue=default         -->" \
          "Wahl der Warteschlange default"
  echo -e "   ABAQUSER.sh terminate             -->" \
          "Vorzeitiges Beenden eines Jobs auf Galilei"
  echo -e "   ABAQUSER.sh suspend               -->" \
          "Unterbrechen eines Jobs auf Galilei"
  echo -e "   ABAQUSER.sh resume                -->" \
          "Wiederaufnehmen eines Jobs auf Galilei"
  echo -e "   ABAQUSER.sh snapshot              -->" \
          "Momentaufnahme eines Jobs auf Galilei"
  echo -e "   ABAQUSER.sh newodb                -->" \
          "Anlegen eines neuen odb nur mit UEL"
  echo -e "   ABAQUSER.sh modifiedupdate        -->" \
          "modifizierte Form (nur in Kombination mit newodb)"
  echo -e "   ABAQUSER.sh urz                   -->" \
          "Ausfuehren des Skriptes auf URZ HPC-Cluster mit qsub"
  echo -e "   ABAQUSER.sh hours=anzahl          -->" \
          "erwartete Dauer der Berechnung in Stunden (nur bei URZ)"
  echo -e "   ABAQUSER.sh mem=speicher          -->" \
          "Hauptspeicher je Chunk (nur bei URZ), z.B. mem=4gb"
  echo -e "   ABAQUSER.sh chunks=anzahl         -->" \
          "Anzahl an Chunks (Nodes) (nur bei URZ), z.B. chunks=2" \
          " Quotient aus Anzahl CPUs und Anzahl Chunks muss ganzzahlig sein !"
  echo -e "Dateiendungen bei input, usr, sharedlib, env, uelinfo mit angeben, Pfadangabe nicht notwendig!"
else
  echo -e "$# Parameter eingegeben"
  echo -e "Eingabedaten:"
  num_opt=0
  k_env=0
  k_qd=0
  k_user=0
  k_sl=0
  k_of=0
  k_cpu=0
  k_job=0
  k_inp=0
  k_info=0
  k_opt=0
  k_fil=0
  k_na=0
  k_replace=0
  k_abaversion=0
  k_mpfil=0
  k_datacheck=0
  k_galilei=0
  k_queue=0
  k_newodb=0
  k_mod=0
  k_terminate=0
  k_suspend=0
  k_resume=0
  k_snapshot=0
  k_oldjob=0
  execstring=""
  usrstring=""
  cpus=1
  k_urz=0
  k_hours=0
  hours=0
  k_mem=0
  memstring=""
  k_chunks=0
  chunks=1
  declare -a koptarray
  declare -a voptarray
  for i in $*
    do
      echo -e "Input: $i"
      k_input=`echo  $i | awk -v FS="=" '{print $1}'`
      v1_input=`echo  $i | awk -v FS="=" '{print $2}'`
      v2_input=`echo  $i | awk -v FS="=" '{print $3}'`
      if [ "$k_input" = "job" ]
        then job="$v1_input"
        jobstring="$v1_input"
        k_job=1
      elif [ "$k_input" = "input" ]
        then
        inputstring="$v1_input"
        k_inp=1
      else
        execstring="$execstring"" ""$i"
        ## Abwaertkompatibilitaet: Optionen uel und umat
        if [ "$k_input" = "usr" ] || [ "$k_input" = "uel" ] || [ "$k_input" = "umat" ]
          then k_user=1
          if [ "$usrstring" = "" ]
            then usrstring="$v1_input"
          else
            usrstring="$usrstring"",""$v1_input"
          fi
        elif [ "$k_input" = "sharedlib" ]
          then slstring="$v1_input"
          k_sl=1
          slAbaqusStandard="libstandardU.so"
        elif [ "$k_input" = "objectfile" ]
          then ofstring="$v1_input"
          k_of=1
          userObjectFile="UserObjectFile.o"
        elif [ "$k_input" = "uelinfo" ]
          then infostring="$v1_input"
          k_info=1
        elif [ "$k_input" = "env" ]
          then envstring="$v1_input"
          k_env=1
        elif [ "$k_input" = "qd" ]
          then
          k_qd=1
        elif [ "$k_input" = "cpus" ]
          then cpus="$v1_input"
          k_cpu=1
        elif [ "$k_input" = "urz" ]
          then k_urz=1
        elif [ "$k_input" = "hours" ]
          then hours="$v1_input"
          k_hours=1
        elif [ "$k_input" = "mem" ]
          then memstring="$v1_input"
          k_mem=1
        elif [ "$k_input" = "chunks" ]
          then chunks="$v1_input"
          k_chunks=1
        elif [ "$k_input" = "opt" ]
          then num_opt=`expr $num_opt + 1`
          k_opt=1
          k_option="$v1_input"
          v_option="$v2_input"
          if [ "$k_option" = "oldjob" ]
           then
            k_oldjob=1
            oldjob="$v_option"
          fi
          koptarray[$num_opt]="$v1_input"
          voptarray[$num_opt]="$v2_input"
        elif [ "$k_input" = "fil" ]
          then k_fil=1
        elif [ "$k_input" = "noanalysis" ]
          then k_na=1
        elif [ "$k_input" = "replace" ]
          then k_replace=1
        elif [ "$k_input" = "abacall" ]
          then abaversion="$v1_input"
          k_abaversion=1
        elif [ "$k_input" = "datacheck" ]
          then k_datacheck=1
        elif [ "$k_input" = "mpfil" ]
          then k_mpfil=1
        elif [ "$k_input" = "galilei" ]
          then k_galilei=1
        elif [ "$k_input" = "queue" ]
          then k_queue=1
          queuestring="$v1_input"
        elif [ "$k_input" = "terminate" ]
          then k_terminate=1
        elif [ "$k_input" = "suspend" ]
          then k_suspend=1
        elif [ "$k_input" = "resume" ]
          then k_resume=1
        elif [ "$k_input" = "snapshot" ]
          then k_snapshot=1
        elif [ "$k_input" = "newodb" ]
          then k_newodb=1
        elif [ "$k_input" = "modifiedupdate" ]
          then k_mod=1
        fi
      fi
  done
  if [ "$k_qd" = 0 ]
  then
    numjobs=0
    if [ "$k_job" = 1 ]
      then numjobs=`number $jobstring`
    fi
    numinps=0
    if [ "$k_inp" = 1 ]
      then numinps=`number $inputstring`
    fi
    numQ=`max $numjobs $numinps`
    for (( i = 1; i <= $numQ; i++ ))
      do
      execstring_i=""
      if [ $i -le $numinps ]
        then
        inp_pfi=`getpathandfile $inputstring $i pf`
        inp_fi=`getpathandfile $inputstring $i f`
        if [ "$inp_fi" != "" ]
          then
          execstring_i="$execstring_i"" ""input=""$inp_pfi"
          name_i=`echo  $(getpathandfile $inputstring $i f) | awk -v FS=".inp" '{print $1}' `
        fi
      fi
      if [ $i -le $numjobs ]
        then
        job_i=`getpathandfile $jobstring $i f`
      elif [ "$name_i" != "" ]
        then
        job_i=$name_i
      else
        job_i="UNNAMEDJOB_""$i"
      fi
      if [ "$job_i" != "" ]
        then
        execstring_i="$execstring_i"" ""job=""$job_i"
        name_i="$job_i"
      fi
      execstring_i="$execstring_i"" ""$execstring"" ""qd"
      if [ \( \( "$k_galilei" -eq 1 \) -o \( "$k_urz" -eq 1 \) \) -a \( "$k_terminate" -eq 0 \) -a \( "$k_suspend" -eq 0 \) -a \( "$k_resume" -eq 0 \) -a \( "$k_snapshot" -eq 0 \) ]
        then
        Qcpus=$cpus
        if [ "$k_na" = 0 ] && [ "$k_mpfil" = 1 ] && [ $cpus < 2]
        then
          Qcpus=2
        elif [ "$k_na" = 1 ] && [ "$k_mpfil" = 0 ]
        then
          Qcpus=1
        elif [ "$k_na" = 1 ] && [ "$k_mpfil" = 1 ]
        then
          Qcpus=2
        fi
        if [ "$k_datacheck" = 1 ]
          then
          Qcpus=1
        fi
        if [ "$k_galilei" = 1 ]
        then
          if [ "$k_queue" = 0 ]
          then
            queue=$queuedefaults
          else
            queue=$queuestring
          fi
          if [ "$k_mem" = 0 ]
          then
            qsub -N "$name_i"_ABAQUSER -j y -q "$queue" -pe shared $Qcpus ABAQUSER.sh $execstring_i
          else
            qsub -N "$name_i"_ABAQUSER -j y -q "$queue" -pe shared $Qcpus -l mem_free="$memstring" ABAQUSER.sh $execstring_i
          fi
        elif [ "$k_urz" = 1 ]
        then
          if [ "$k_queue" = 0 ]
          then
            queue=""
          else
            queue="-q $queuestring"
          fi
          if [ "$k_datacheck" = 1 ]
            then
            chunks=1
          fi
          batchJob="ABAQUSER-Job"
          NcpusPerChunk=`expr $Qcpus / $chunks`
          qsubOpts="-l select="$chunks":ncpus="$NcpusPerChunk":mem="$memstring":mpiprocs="$NcpusPerChunk" -l walltime="$hours":00:00 -o "$job".log -j oe "$queue" -N $batchJob"
          echo "$qsubOpts"
          SCRIPTNAME="ABAQUSER.sh"
          echo "$PWD/$SCRIPTNAME $execstring_i" | qsub $qsubOpts
        fi
      else
        ./ABAQUSER.sh $execstring_i
      fi
    done
  else
    if [ "$k_abaversion" = 1 ]
    then
      aba=$abaversion
    else
      aba=$latestabacall
    fi
    if [ "$k_urz" = 1 ]
    then
      echo "Loading modules"
      echo "change abacall to abaqus"
      aba="abaqus"
      if [ "$abaversion" = "abq610ef1" ]
        then
        module load abaqus/6.10
      elif [ "$abaversion" = "abq6123" ]
        then
        module load abaqus/6.12
      elif [ "$abaversion" = "abq6142" ]
        then
        module load abaqus/6.14
     else
        echo "No ABAQUS module available. Choose other abacall:" $abaversion
      fi
    fi
    if [ "$k_inp" = 1 ]
    then
      inputname=`echo  $(getpathandfile $inputstring 1 f) | awk -v FS=".inp" '{print $1}' `
      cp $job.inp $job"_old".inp_old 2>> out_$job.dat
      cp $inputstring $job.inp 2>> out_$job.dat
      if [ "$k_na" = 1 ]
      then
        cp $inputname.fil $job.fil 2>> out_$job.dat
        cp $inputname.odb $job.odb 2>> out_$job.dat
      fi
    fi
    if [ "$k_replace" = 1 ] && [ "$k_na" = 0 ] && [ "$k_galilei" = 0 ]
      then
      mv $job.inp inputtemp_$job.inp
      rm $job.* 2>> out_$job.dat
      mv inputtemp_$job.inp $job.inp
    fi
    if [ "$k_na" = 1 ]
      then
      cp $job.odb $job"_SDV".odb
      cp $job.fil $job"_SDV".fil
      cp $job.inp $job"_SDV".inp
      jobold=$job
      job=$job"_SDV"
    fi
    if [ "$k_na" = 0 ]
      then
      if [ "$k_fil" = 1 ]
        then
        if [ "$k_urz" = 0 ]
          then
          /app/bin/$aba python -c "__import__('ABAQUSER').InputFile('$job')"
        fi
      fi
      if [ "$k_user" = 1 ] && [ "$k_sl" = 1 ]
        then
        echo -e "Shared library is specified! Option usr is ignored."
      elif [ "$k_user" = 1 ] && [ "$k_of" = 1 ]
        then
        echo -e "Object file is specified! Option usr is ignored."
      elif [ "$k_user" = 1 ] && [ "$k_sl" = 0 ] && [ "$k_of" = 0 ]
        then
        rm userfile_$job.f 2>> out_$job.dat
        numusr=`number $usrstring`
        for (( filenum = 1; filenum <= $numusr; filenum++ ))
          do
          pathfile=`getpathandfile $usrstring $filenum pf`
          cat "$pathfile" >> userfile_$job.f
        done
      fi
      if [ "$k_job" = 0 ]
        then
        echo -e "Kein Job definiert! " \
                "Erneute Eingabe mit abaquser job=jobname erforderlich!"
      else
        if [ "$k_datacheck" = 1 ]
          then
          echo -e "DATACHECK"
          abacall="$aba"" job=""$job"" ""datacheck"" ""interactiv"
        else
          abacall="$aba"" job=""$job"" ""analysis"" ""interactiv"
        fi
        if [ "$k_of" = 1 ] && [ "$k_sl" = 0 ]
          then   abacall="$abacall"" ""user=$userObjectFile"
        elif [ "$k_user" = 1 ] && [ "$k_sl" = 0 ]
          then   abacall="$abacall"" ""user=userfile_$job.f"
        fi
        if [ "$k_datacheck" = 1 ]
          then
          cpus=1
        fi
        if [ "$k_cpu" = 1 ]
          then   abacall="$abacall"" ""cpus=""$cpus"
        fi
        if [ "$k_opt" = 1 ]
          then
          for (( i = 1; i <= $num_opt; i++ ))
            do
            abacall="$abacall"" ""${koptarray[$i]}""=""${voptarray[$i]}"
          done
        fi
        echo -e "Aufruf von ABAQUS mit:"
        echo -e "$abacall"
        date
        if [ "$k_env" = 1 ]
          then
          envpathfile=`getpathandfile $envstring 1 pf`
          envfile=`getpathandfile $envstring 1 f`
        fi
        if [ "$k_sl" = 1 ]
          then
          slpathfile=`getpathandfile $slstring 1 pf`
          slfile=`getpathandfile $slstring 1 f`
        fi
        if [ "$k_of" = 1 ]
          then
          ofpathfile=`getpathandfile $ofstring 1 pf`
          offile=`getpathandfile $ofstring 1 f`
        fi
        if [ "$k_galilei" = 0 ] && [ "$k_urz" = 0 ]
          then
          if [ "$k_env" = 1 ]
            then
            cp "$envpathfile" ./
          fi
          if [ "$k_sl" = 1 ]
            then
            cp "$slpathfile" ./$slAbaqusStandard
          fi
          if [ "$k_of" = 1 ]
            then
            cp "$ofpathfile" ./$userObjectFile
          fi
          $abacall
        elif [ "$k_galilei" = 1 ] && [ "$k_terminate" = 0 ] && [ "$k_suspend" = 0 ] && [ "$k_resume" = 0 ] && [ "$k_snapshot" = 0 ]
          then
          #$ -pe shared $cpus
          echo -e "Anzahl CPUS:" $cpus
          EXEC=$abacall
          #$ -cwd
          FILES=`ls -1 ${job}.inp *.f *.py *.pyc *.aba *.F *.f90 *.o 2>> out_$job.dat`
          if [ "$k_env" = 0 ]
            then
            FILES=`ls -1  $FILES *.env 2>> out_$job.dat`
          fi
          if [ "$k_oldjob" = 1 ]
            then
            FILES=`ls -1  $FILES ${oldjob}.* 2>> out_$job.dat`
          fi
          EXEC="$EXEC"
          NODE=`hostname`
          SDIR=$PWD
          NDIR=/scratch/$USER.$JOB_ID
          export PATH=$PATH:/app/bin
          export TMPDIR=$NDIR
          mkdir $NDIR
          chmod 750 $NDIR
          ln -s /scratch/$NODE/$USER.$JOB_ID  $SDIR/scratch.$JOB_ID
          if [ -n "$FILES" ] ; then
            cp $FILES $NDIR
            if [ "$k_env" = 1 ]
              then
              cp $envpathfile $NDIR
            fi
            if [ "$k_sl" = 1 ]
              then
              cp $slpathfile $NDIR/$slAbaqusStandard
            fi
            if [ "$k_of" = 1 ]
              then
              cp $ofpathfile $NDIR/$userObjectFile
            fi
          fi
          cd $NDIR
          echo "Start checkFCGProgress"
          checkFCGProgress_execstring="__import__('FCGProgress').checkFCGProgress(['job=$job', 'abacall=$aba', 'sleeptime=3600'])"
          echo $checkFCGProgress_execstring
          $aba python -c "$checkFCGProgress_execstring" 1>$job.check 2>$job.errcheck &
          $EXEC
          if [ -n "$FILES" ] ; then
            rm $FILES 2>> out_$job.dat
          fi
          if [ "$k_env" = 1 ]
            then
            rm $envfile 2>> out_$job.dat
          fi
          if [ "$k_sl" = 1 ]
            then
            rm $slAbaqusStandard 2>> out_$job.dat
          fi
          if [ "$k_of" = 1 ]
            then
            rm $userObjectFile 2>> out_$job.dat
          fi
          if [ "$k_fil" = 0 ] && [ "$k_mpfil" = 0 ]
            then
            cp ./* $SDIR
            cd ..
            rm -r ./$USER.$JOB_ID
            rm $SDIR/scratch.$JOB_ID
            echo Job $JOB_ID beendet
          fi
          cd $SDIR
        elif [ "$k_urz" = 1 ] && [ "$k_terminate" = 0 ] && [ "$k_suspend" = 0 ] && [ "$k_resume" = 0 ] && [ "$k_snapshot" = 0 ]
          then
          echo -e "Anzahl CPUS:" $cpus
          EXEC=$abacall
          WORKDIR="$PANFS/$PBS_JOBID"
          LINKNAME="panfs_$PBS_JOBID"
          NCPUSTOT=`qstat -f $PBS_JOBID | sed -n -e 's/ //g' -e 's/Resource_List.ncpus=//p'`
          mkdir "$WORKDIR"
          cd "$PBS_O_WORKDIR"
          ln -s "$WORKDIR" "$LINKNAME"
          FILES=`ls -1 ${job}.inp *.f *.py *.pyc *.aba *.F *.f90 *.o 2>> out_$job.dat`
          if [ "$k_env" = 0 ]
            then
            FILES=`ls -1  $FILES *.env 2>> out_$job.dat`
          fi
          if [ "$k_oldjob" = 1 ]
            then
            FILES=`ls -1  $FILES ${oldjob}.* 2>> out_$job.dat`
          fi
          EXEC="$EXEC"
          if [ -n "$FILES" ] ; then
            cp $FILES $WORKDIR
            if [ "$k_env" = 1 ]
              then
              cp $envpathfile $WORKDIR
            fi
            if [ "$k_sl" = 1 ]
              then
              cp $slpathfile $WORKDIR/$slAbaqusStandard
            fi
            if [ "$k_of" = 1 ]
              then
              cp $ofpathfile $WORKDIR/$userObjectFile
            fi
          fi
          cd $WORKDIR
          if [ "$k_fil" = 1 ]
            then
            $aba python -c "__import__('ABAQUSER').InputFile('$job')"
          fi
          echo "Start checkFCGProgress"
          checkFCGProgress_execstring="__import__('FCGProgress').checkFCGProgress(['job=$job', 'abacall=$aba', 'sleeptime=3600'])"
          echo $checkFCGProgress_execstring
          $aba python -c "$checkFCGProgress_execstring" 1>$job.check 2>$job.errcheck &
          echo "Start ABAQUS"
          echo -e $EXEC
          $EXEC
          rm *.cluster 2>> out_$job.dat
          if [ -n "$FILES" ] ; then
            rm $FILES 2>> out_$job.dat
          fi
          if [ "$k_env" = 1 ]
            then
            rm $envfile 2>> out_$job.dat
          fi
          if [ "$k_sl" = 1 ]
            then
            rm $slAbaqusStandard 2>> out_$job.dat
          fi
          if [ "$k_of" = 1 ]
            then
            rm $userObjectFile 2>> out_$job.dat
          fi
          if [ "$k_fil" = 0 ] && [ "$k_mpfil" = 0 ]
            then
            mv ./* $PBS_O_WORKDIR 2>> out_$job.dat
            cd "$PBS_O_WORKDIR"
            rmdir "$WORKDIR" 2>> out_$job.dat
            rm "$LINKNAME" 2>> out_$job.dat
            echo Job $JOB_ID beendet
          fi
          cd $PBS_O_WORKDIR
        fi
        date
      fi
      if [ "$k_user" = 1 ] && [ "$k_sl" = 0 ] && [ "$k_of" = 0 ]
        then rm userfile_$job.f 2>> out_$job.dat
      fi
    fi
    if [ \( \( "$k_fil" -eq 1 \) -o \( "$k_mpfil" -eq 1 \) \) -a \( "$k_terminate" -eq 0 \) -a \( "$k_suspend" -eq 0 \) -a \( "$k_resume" -eq 0 \) -a \( "$k_snapshot" -eq 0 \) ]
      then
      if [ "$k_mpfil" = 1 ] && [ $aba != $latestabacall ]
      then
        if [ $k_replace = 1 ]
        then
          rm $job"_"$latestabacall.* 2>> out_$job.dat
        fi
        echo -e "START ODB FILE UPGRADE" $aba "-->" $latestabacall
        "/app/bin/"$latestabacall upgrade job=$job"_"$latestabacall odb=$job
        #$upgradestring
        cp $job.fil $job"_"$latestabacall.fil 2>> out_$job.dat
        jobold2=$job
        job="$job""_""$latestabacall"
        aba=$latestabacall
        rm "$job"-upgrade.log 2>> out_$job.dat
      fi
      fil2odb_execstring="['$job'"",'ABACALL=""$aba""'"
      if [ "$k_newodb" = 1 ]
        then fil2odb_execstring="$fil2odb_execstring"",'NEWODB'"
      fi
      if [ "$k_newodb" = 1 ] && [ "$k_mod" = 1 ]
        then fil2odb_execstring="$fil2odb_execstring"",'MODIFIEDUPDATE'"
      fi
      if [ "$k_mpfil" = 1 ]
        then fil2odb_execstring="$fil2odb_execstring"",'MP'"
      fi
      if [ "$k_info" = 0 ]
        then infopathfile="./UEL.info"
        infofile="UEL.info"
        echo -e "Kein *.info definiert. Versuch mit UEL.info..."
      else
        infopathfile=`getpathandfile $infostring 1 pf`
        infofile=`getpathandfile $infostring 1 f`
      fi
      if [ "$k_galilei" = 0 ] && [ "$k_urz" = 0 ]
        then
        fil2odb_execstring="$fil2odb_execstring"",'$infopathfile'"
        fil2odb_execstring="__import__('ABAQUSER').Fil2Odb($fil2odb_execstring])"
        echo -e "Aufruf von fil2odb mit:"
        echo -e "$fil2odb_execstring"
        $aba python -c "$fil2odb_execstring"
      elif [ "$k_galilei" = 1 ]
        then
        cp "$infopathfile" ./
        fil2odb_execstring="$fil2odb_execstring"",'$infofile'"
        fil2odb_execstring="__import__('ABAQUSER').Fil2Odb($fil2odb_execstring])"
        if [ "$k_mpfil" = 0 ]
          then
          numcpu=1
        elif [ "$k_mpfil" = 1 ]
          then
          numcpu=2
        fi
        #$ -pe shared $numcpu
        echo -e "Anzahl CPUS:" $numcpu
        EXEC="/app/bin/$aba python -c "'"'$fil2odb_execstring'"'
        echo -e "Aufruf von fil2odb:"
        echo -e $EXEC
        #$ -cwd
        if [ "$k_na" = 1 ]
        then
          FILES1=`ls -1 ${job}.odb 2>> out_$job.dat`
          FILES2=`ls -1 ${job}.fil *.py* $infofile 2>> out_$job.dat`
        else
          FILES2=`ls -1 *.py* $infofile 2>> out_$job.dat`
        fi
        EXEC="$EXEC"
        NODE=`hostname`
        SDIR=$PWD
        NDIR=/scratch/$USER.$JOB_ID
        export PATH=$PATH:/app/bin
        export TMPDIR=$NDIR
        mkdir $NDIR
        chmod 750 $NDIR
        ln -s /scratch/$NODE/$USER.$JOB_ID  $SDIR/scratch.$JOB_ID
        if [ "$k_na" = 1 ] && [ -n "$FILES1" ] ; then
          cp $FILES1 $NDIR
        fi
        if [ -n "$FILES2" ] ; then
          cp $FILES2   $NDIR
        fi
        cd $NDIR
        eval $EXEC
        if [ -n "$FILES2" ] ; then
          rm $FILES2 2>> out_$job.dat
        fi
        if [ "$k_na" = 0 ]
        then
          rm ${job}.fil 2>> out_$job.dat
          rm *.pyc 2>> out_$job.dat
        fi
        if [ "$k_na" = 1 ]
        then
          cp $FILES1 $SDIR 2>> out_$job.dat
        else
          cp ./* $SDIR 2>> out_$job.dat
        fi
        cd ..
        rm -r ./$USER.$JOB_ID
        rm $SDIR/scratch.$JOB_ID
        echo Job $JOB_ID beendet
        cd $SDIR
      elif [ "$k_urz" = 1 ]
        then
        fil2odb_execstring="$fil2odb_execstring"",'$infofile'"
        fil2odb_execstring="__import__('ABAQUSER').Fil2Odb($fil2odb_execstring])"
        if [ "$k_mpfil" = 0 ]
          then
          numcpu=1
        elif [ "$k_mpfil" = 1 ]
          then
          numcpu=2
        fi
        echo -e "Anzahl CPUS:" $numcpu
        EXEC="$aba python -c "'"'$fil2odb_execstring'"'
        echo -e "Aufruf von fil2odb:"
        echo -e $EXEC
        WORKDIR="$PANFS/$PBS_JOBID"
        LINKNAME="panfs_$PBS_JOBID"
        mkdir "$WORKDIR"
        cd "$PBS_O_WORKDIR"
        ln -s "$WORKDIR" "$LINKNAME"
        cp "$infopathfile" $WORKDIR
        if [ "$k_na" = 1 ]
        then
          FILES1=`ls -1 ${job}.odb 2>> out_$job.dat`
          FILES2=`ls -1 ${job}.fil *.py* $infofile 2>> out_$job.dat`
        else
          FILES2=`ls -1 *.py* $infofile 2>> out_$job.dat`
        fi
        EXEC="$EXEC"
        if [ "$k_na" = 1 ] && [ -n "$FILES1" ] ; then
          cp $FILES1 $WORKDIR
        fi
        if [ -n "$FILES2" ] ; then
          cp $FILES2 $WORKDIR
        fi
        cd $WORKDIR
        eval $EXEC
        if [ -n "$FILES2" ] ; then
          rm $FILES2 2>> out_$job.dat
        fi
        if [ "$k_na" = 0 ]
        then
          rm ${job}.fil 2>> out_$job.dat
          rm *.pyc 2>> out_$job.dat
        fi
        rm *.cluster 2>> out_$job.dat
        if [ "$k_na" = 1 ]
        then
          mv $FILES1 $PBS_O_WORKDIR 2>> out_$job.dat
        else
          mv ./* $PBS_O_WORKDIR 2>> out_$job.dat
        fi
        cd "$PBS_O_WORKDIR"
        rmdir "$WORKDIR" 2>> out_$job.dat
        rm "$LINKNAME" 2>> out_$job.dat
        echo Job $JOB_ID beendet
        cd $PBS_O_WORKDIR
      fi
      date
    fi
    if [ \( "$k_galilei" -eq 1 \) -a \( \( "$k_terminate" -eq 1 \) -o \( "$k_suspend" -eq 1 \) -o \( "$k_resume" -eq 1 \) -o \( "$k_snapshot" -eq 1 \) \) ]
      then
      job2kill="$job"
      jobnameinqueue="$job2kill""_ABAQUSER"
      jobtest=`qstat`
      if [ "$jobtest" != "" ]
        then
        jobids=`qstat -j $jobnameinqueue | sed -n -e  's/job_number://p'`
        jobnames=`qstat -j $jobnameinqueue | sed -n -e  's/job_name://p'`
        jobargs=`qstat -j $jobnameinqueue | sed -n -e  's/job_args://p'`
        jobid=`echo  $jobids | awk -F'\n' '{print $0}'`
        jobname=`echo  $jobnames | awk -F'\n' '{print $0}'`
        jobarg=`echo  $jobargs | awk -F'\n' '{print $0}'`
        numjob=`echo $jobid | awk -F' ' '{print NF}'`
        k_job=0
        k_kill=0
        id_counter=0
        while [ "$k_job" = 0 ]
          do
          optvalue="notthe""$job2kill"
          id_counter=$[$id_counter+1]
          id=`echo $jobid | awk -F' ' '{print $id_counter}' id_counter=$id_counter`
          name=`echo $jobname | awk -F' ' '{print $id_counter}' id_counter=$id_counter`
          arg=`echo $jobarg | awk -F' ' '{print $id_counter}' id_counter=$id_counter`
          if [ "$name" = "$jobnameinqueue" ]
            then
            k_jobopt=0
            counter=0
            while [ "$k_jobopt" = 0 ]
              do
              counter=$[$counter+1]
              opt=`echo  $arg | awk -F',' '{print $counter}' counter=$counter`
              optname=`echo $opt | awk -F'=' '{print $1}'`
              optvalue=`echo $opt | awk -F'=' '{print $2}'`
              if [ "$optname" = "job" ]
                then
                k_jobopt=1
              fi
            done
          elif [ "$id_counter" = $numjob ]
            then
            echo "No job" "$jobnameinqueue" "found in queue!"
            k_job=1
          fi
          if [ "$optvalue" = "$job2kill" ]
          then
            k_job=1
            id2kill=$id
            k_kill=1
          elif [ "$name" = "$jobnameinqueue" ] && [ "$id_counter" = $numjob ]
          then
            echo "Job" "$job2kill" "not found!"
            k_job=1
          fi
        done
        if [ "$k_kill" = 1 ]
        then
          node_ext=`qstat | sed -n -e  s/$id2kill//p | sed -n -e s/.\*@//p`
          node=`echo $node_ext | awk -F' ' '{print $1}'`
          if [ "$k_terminate" = 1 ]
          then
            command="cd /scratch/$USER.$id2kill/;rm $job2kill.lck;/app/bin/$aba job=$job2kill suspend;sleep 10;/app/bin/$aba job=$job2kill terminate;exit"
            ssh $node $command
            sleep 5
          elif [ "$k_suspend" = 1 ]
          then
            command="cd /scratch/$USER.$id2kill/;rm $job2kill.lck;/app/bin/$aba job=$job2kill suspend;exit"
            ssh $node $command
            sleep 5
          elif [ "$k_resume" = 1 ]
          then
            command="cd /scratch/$USER.$id2kill/;rm $job2kill.lck;/app/bin/$aba job=$job2kill resume;exit"
            ssh $node $command
            sleep 5
          elif [ "$k_snapshot" = 1 ]
          then
            command="cd /scratch/$USER.$id2kill/;/app/bin/$aba job=$job2kill suspend;exit"
            ssh $node $command
            sleep 5
            timestamp=`date +"%y%m%d%H%M%S"`
            command="cp /scratch/$USER.$id2kill/$job2kill.odb $PWD/$job2kill-$timestamp.odb"
            ssh $node $command
            sleep 5
            command="cp /scratch/$USER.$id2kill/$job2kill.fil $PWD/$job2kill-$timestamp.fil"
            ssh $node $command
            sleep 5
            echo snapshot $job2kill-$timestamp created
            command="cd /scratch/$USER.$id2kill/;/app/bin/$aba job=$job2kill resume;exit"
            ssh $node $command
            sleep 5
          fi
        fi
      else
        echo "No jobs found in queue!"
      fi
    fi
    rm -f $job.{023,abq,com,ipm,pac,par,pes,pmg,sel,lck,023,cid,odb_f} 2>> out_$job.dat
    if [ "$k_na" = 1 ]
      then
      rm -f $job.{fil,inp} 2>> out_$job.dat
    fi
    rm -f out_$jobold{.dat,2.dat} 2>> out_$job.dat
    rm out_$job.dat
    echo -e "\nABAQUSER SUCCESSFULLY COMPLETED"
  fi
fi
##########################################################################
