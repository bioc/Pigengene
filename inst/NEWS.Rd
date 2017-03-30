\name{NEWS}
\title{News for \R Package \pkg{Pigenegene}}
   
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
