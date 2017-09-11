\name{menina-package}
\alias{menina-package}
\alias{menina}
\docType{package}
\title{
  menina implementa a menor infraestrutura de acesso em redes sem fio obstrui'das.
}
\description{
  O pacote permite:
  Definir a instancia do problema: um grafo imerso em um grid bidimensional.
  Calcular os componentes do grafo da rede sem fio obstruida.
  Buscar os pontos base candidatos a colocacao de pontos de acesso.
  Assentar nos melhores pontos base os acessos pelo uso do algoritmo menina e do guloso.
}
\details{
  D = definaDispositivo(g, mu, modo)
  C = computeComponente(D, r, eps)
  B = busqueBase(D, C, g, r, eps)
  A = assenteAcesso(C, B, algoritmo)
}
\author{
Manasses Ferreira

Maintainer: Manasses Ferreira <manassesferreira@gmail.com>
}
\references{
  AUTHOR="Manasses Ferreira, Neto and Olga Goussevskaia and Vinicius {dos Santos}",
  TITLE="Connectivity with Backbone Structures in Obstructed Wireless Networks",
  BOOKTITLE="COMPUTER NETWORKS",
  YEAR=2017
}
\keyword{ package, analytical modeling and model validation, performance,
  fundamental limits, scalability, energy, and reliability, communication networks
 }
\seealso{
  Optional links to other man pages
}
\examples{
  \dontrun{
     ## Optional simple examples of the most important functions
     ## These can be in \dontrun{} and \donttest{} blocks.
  }
}
