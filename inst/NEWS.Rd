\name{NEWS}
\title{News for \R Package \pkg{Pigenegene}}
   
\section{Changes in version 1.3.6 (2017-08-20)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item A bug in gene.mapping () function fixed to better map probe IDs.
    }
  }
}

\section{Changes in version 1.3.4 (2017-07-31)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The doTranspose argument added to the  heatmap.type() function.
    }
  }
}

\section{Changes in version 1.1.14 (2017-05-10)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item A bug in the module.heatmap() function fixed.
    }
  }
}

\section{Changes in version 1.1.12 (2017-04-15)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The mouseHomologFilter argument is added to the unexported
      furntion hu.mouse() with a default value of
      "with_mmusculus_homolog". This is now compatible with
      biomaRt_2.30.0, it used to be "with_homolog_mmus" in older versions of biomaRt.
    }
  }
}

\section{Changes in version 1.1.10 (2017-03-30)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item Checking the pigengene input of module.heatmap().
      \item Issues in the balance() function (not exported) where Labels
      is a factor are resolved. Also,  if all sampls have the same size,
      oversampling is automatically turned off.
      \item If Labels is a factor, it is now converted to a character
      vector in check. pigengene.input().
    }
  }
}

\section{Changes in version 1.1.6 (2017-03-27)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The module.heatmap() function now has the doAddEigengene
      and scalePngs arguments.
      \item The compute.pigengene()  function now reports also the size of
      modules in the pigengene_pvalue.csv output file.
    }
  }
}

\section{Changes in version 0.99.25 (2016-10-02)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The compute.pigengene() function now uses welch.pvalue() instead of pvalues.manov().
    }
  }
}

\section{Changes in version 0.99.8 (2016-05-18, git f50187c)}{
  \subsection{General}{
    \itemize{
      \item  Under review by Bioconductor.
      \item Created.
    }
  }
}

%% \subsection{General}
%% \subsection{New functions}{
%% \subsection{Changes in existing functions} %% functionality, signature, output
%% \subsection{Bug Fixes}
