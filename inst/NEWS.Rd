\name{NEWS}
\title{News for \R Package \pkg{Pigenegene}}


\section{Changes in version 1.11.34 (2019-10-21)}{
  \subsection{General}{
    \itemize{
      \item The pipeline is now explained step by step in the vignette.  
    }
  }
}

\section{Changes in version 1.11.32 (2019-10-18)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The \code{Data} argument of \code{compute.pigengene} can now
      be a matrix with only 1 column. 
    }
  }
}

\section{Changes in version 1.11.30 (2019-10-02)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The \code{doRetuNetworks} argument is now added to the
      \code{combine.networks} function.
    }
  }
}

\section{Changes in version 1.11.28 (2019-10-01)}{
  \subsection{New functions}{
    \itemize{
      \item \code{message.if()} is now exported.
    }
  }
}

\section{Changes in version 1.11.26 (2019-09-26)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item Better QC in the \code{gene.mapping()} function, the
      possible keys will be printed if the input is not appropriate.
    }
  }
  \subsection{New functions}{
    \itemize{
      \item \code{save.if()} is now exported.
    }
  }
}

\section{Changes in version 1.11.24 (2019-09-03)}{
  \subsection{Bug Fixes}{
    \itemize{
      \item In the \code{combine.network()} function,
      \code{selectedModules} does not need to be \code{"All"}. Also, if
      \code{saveFile=NULL}, nothing will be saved without any error.
    }
  }
}

\section{Changes in version 1.11.20 (2019-05-15)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item \code{Data} and \code{Labels} can now be lists, which will be combined using 
      \code{combine.network()} before analysis.
    }
  }
}

\section{Changes in version 1.11.4 (2019-05-02)}{
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{repeat.data(times=1,...)} now produces valid output.
    }
  }
}

\section{Changes in version 1.9.26 (2019-04-24)}{
  \subsection{Bug Fixes}{
    \itemize{
      \item In the \code{compute.pigengene} function, if a module has
      only one gene, its name is now not omitted in the csv file.
    }
  }
}

\section{Changes in version 1.9.24 (2019-04-24)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item An object of \code{\link{pigengene-class}} can now have a
      \code{heavyToLow} attribute.
    }
  }
}

\section{Changes in version 1.9.20 (2019-04-12)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item In the \code{compute.pigengene} function, the columns of
      \code{Data} and names of \code{modules} can now differ.
    }
  }
}

\section{Changes in version 1.9.14 (2019-02-12)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The \code{dOrderByW} argument is now added for the \code{compute.pigengene} function.
    }
  }
}

\section{Changes in version 1.9.8 (2018-11-22)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item All of the \code{learn.bn} function can now be set from
      \code{one.step.pigengene} through \code{bnArgs}.
    }
  }
}

\section{Changes in version 1.9.4 (2018-11-16)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item \code{gene.mapping} can now output multiple conventions.
    }
  }
}

\section{Changes in version 1.7.2 (2018-05-22)}{
  \subsection{General}{
    \itemize{
      \item The version of the package \code{C50} is now required to be
      at least 0.1.2, which exports the \code{as.party.C5.0()} function.
    }
  }
}

\section{Changes in version 1.5.22 (2018-04-28)}{
  \subsection{Bug Fixes}{
   \itemize{
      \item \code{stats::cor} is used in \code{compute.pigengene} and
  \code{draw.cor} functions. See the NAMESPACE for the important reason.
    }
  }
}

\section{Changes in version 1.5.9 (2018-03-12)}{
  \subsection{General}{
    \itemize{
      \item \code{cor} is imported from WGCNA, but not from stats, because WGCNA does
      not call the \code{cor} function properly.
    }
  }
  \subsection{Changes in existing functions}{
    \itemize{
      \item \code{RsquaredCut} is added to the arguments of the \code{one.step.pigengene} function.
    }
  }
}

\section{Changes in version 1.5.6 (2018-01-19)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item In the \code{get.fitted.leaf} function, the function
  C50:::as.party.C5.0 is used, which used to be exported in the previous
  versions of C50, but not in version 0.1.1.
    }
  }
}

\section{Changes in version 1.5.2 (2017-11-10)}{
  \subsection{New functions}{
    \itemize{
      \item The \code{check.nas} function is now exported.
    }
  }
}

\section{Changes in version 1.3.8 (2017-09-13)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item Order of conditions in pheatmap.type can now be determined
  by the user.
    }
  }
}

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
