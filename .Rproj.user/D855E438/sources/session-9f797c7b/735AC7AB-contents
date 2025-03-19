if(!require('rsconnect')) install.packages('rsconnect'); library('rsconnect')

### conectar no servidor
rsconnect::setAccountInfo(
  name='lipaufabc',
  token='DB04322E5C8F58C48C78A2CF346B220C',
  secret='pHMwoVIy7SKkrxYXoKNpntAHNKswtc2wN4Gg8C0A'
  )
### containerizar a imagem
rsconnect::deployApp(
  appDir = 'C:/Users/eduar/OneDrive/Desktop/GitHub/lipa_dashboard/APP',
  appId = '14225081',
  appName = 'lipa_dashboard',
  account = 'lipaufabc'
)

