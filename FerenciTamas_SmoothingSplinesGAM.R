## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(fig.align="center")
# knitr::opts_chunk$set(fig.align="center", out.width = "70%")
# download.file(paste0("https://raw.githubusercontent.com/tamas-ferenci/FerenciTamas_SmoothingSplinesGAM/",
#                      "master/FerenciTamas_SmoothingSplinesGAM_cover.png"),
#               "./docs/FerenciTamas_SmoothingSplinesGAM_cover.png", mode = "wb")


## Téma: simítás, spline-regresszió, additív modellek


## Minden visszajelzést örömmel veszek a [tamas.ferenci@medstat.hu](tamas.ferenci@medstat.hu) email-címen


## A LOESS simítóról lesz szó


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
set.seed(1)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
n <- 101
x <- (1:n) + rnorm(n, 0, 0.1)
y <- sin(x/n*(2*pi))
yobs <- y + rnorm(n, 0, 0.2)
SimData <- data.frame(x, y, yobs)
p <- ggplot(SimData, aes(x = x, y = yobs)) + geom_point() +
  geom_line(aes(y = y), color = "orange", lwd = 1)
p


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
p + geom_smooth(formula = y~x, method = "lm", se = FALSE)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
p + geom_vline(xintercept = 23.5, color = "red")


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
span <- 0.75
n*span
ceiling(n*span)
sort(abs(x-23.5))
sort(abs(x-23.5))[ceiling(n*span)]


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
tricube <- function(u, t) ifelse(u<t, (1-(u/t)^3)^3, 0)
curve(tricube(x, 2), to = 3)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
SimData$w <- tricube(abs(x-23.5), sort(abs(x-23.5))[ceiling(n*span)])
ggplot(SimData, aes(x = x, y = yobs, color = w>0)) + geom_point()


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
ggplot(SimData, aes(x = x, y = yobs, color = w)) + geom_point()


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
fit <- lm(yobs ~ x, weights = w, data = SimData)
p + geom_vline(xintercept = 23.5, color = "red") +
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2])


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
p + geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2]) +
  geom_vline(xintercept = 23.5, color = "red") +
  geom_point(x = 23.5, y = predict(fit, data.frame(x = 23.5)), color="red")


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
loessfun <- function(xin, x, yobs, span) {
  n <- length(x)
  w <- tricube(abs(x-xin), sort(abs(x-xin))[ceiling(n*span)])
  fit <- lm(yobs ~ x, weights = w)
  predict(fit, data.frame(x = xin))
}
p + geom_vline(xintercept = 23.5, color = "red") +
  geom_point(x = 23.5, y = loessfun(23.5, SimData$x, SimData$yobs, 0.75), color="red")


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
p + geom_vline(xintercept = 48.3, color = "red") +
  geom_point(x = 48.3, y = loessfun(48.3, SimData$x, SimData$yobs, 0.75), color="red")
p + geom_vline(xintercept = 91.2, color = "red") +
  geom_point(x = 91.2, y = loessfun(91.2, SimData$x, SimData$yobs, 0.75), color="red")


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
SmoothData <- expand.grid(grid = seq(0, 101, 0.1))
SmoothData$value <- apply(SmoothData, 1, function(x)
  loessfun(x["grid"], SimData$x, SimData$yobs, 0.75))
p + geom_line(data = SmoothData, aes(x = grid, y = value), color = "red")


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
SmoothData <- expand.grid(grid = seq(0, 101, 0.1),
                          span = c(2/n+1e-10, 0.25, 0.5, 0.75, 1))
SmoothData$value <- apply(SmoothData, 1, function(x)
  loessfun(x["grid"], SimData$x, SimData$yobs, x["span"]))
SmoothData$span <- as.factor(SmoothData$span)

p + geom_line(data = SmoothData[SmoothData$span%in%c(0.25, 0.5, 0.75), ],
              aes(x = grid, y = value, color = span))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
p + geom_line(data = SmoothData[SmoothData$span==1, ], aes(x = grid, y = value),
              color = "red")


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
p + geom_line(data = SmoothData[SmoothData$span==2/n+1e-10, ], aes(x = grid, y = value),
              color = "red")


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
lm(y~x+x^2, data = SimData)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
lm(y~x+I(x^2), data = SimData)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
lm(y~poly(x,2), data = SimData)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
predict(lm(y~x+I(x^2)), data.frame(x = 43.9))
predict(lm(y~poly(x,2)), data.frame(x = 43.9))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
t(cbind(1, poly(x, 3)))%*%cbind(1, poly(x, 3))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
lm(y~poly(x,2, raw = TRUE), data = SimData)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
kappa(cbind(1, poly(x, 3)), exact = TRUE)
kappa(cbind(1, poly(x, 3, raw = TRUE)), exact = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
loessfun <- function(xin, x, yobs, span, degree) {
  n <- length(x)
  w <- tricube(abs(x-xin), sort(abs(x-xin))[ceiling(n*span)])
  fit <- lm(yobs ~ poly(x, degree), weights = w)
  predict(fit, data.frame(x = xin))
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
SmoothData <- rbind(expand.grid(grid = seq(0, 101, 0.1),
                                span = c(2/n+1e-10, 0.25, 0.5, 0.75, 1), degree = 1),
                    expand.grid(grid = seq(0, 101, 0.1),
                                span = c(3/n+1e-10, 0.25, 0.5, 0.75, 1), degree = 2))
SmoothData$value <- apply(SmoothData, 1, function(x)
  loessfun(x["grid"], SimData$x, SimData$yobs, x["span"], x["degree"]))
SmoothData$span <- as.factor(SmoothData$span)
p + geom_line(data = SmoothData[SmoothData$span%in%c(0.25, 0.5, 0.75), ],
              aes(x = grid, y = value, color = span)) + facet_grid(rows = vars(degree))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
SimData$train <- FALSE
SimData$train[sample(1:101, 80)] <- TRUE
ggplot(SimData, aes(x = x, y = yobs, color = train)) + geom_point() +
  geom_line(aes(y = y), color = "orange", lwd = 1)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
spans <- seq(0.03, 0.99, 0.01)
SSEfull <- sapply(spans, function(sp) sum((sapply(1:101, function(i)
  loessfun(SimData$x[i], SimData$x, SimData$yobs, sp, 1))-yobs)^2))
ggplot(data.frame(span = spans, SSE = SSEfull), aes(x = span, y = SSE)) + geom_line()


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
SSEfull <- sapply(spans, function(sp) sum((sapply(which(!SimData$train), function(i)
  loessfun(SimData$x[i], SimData$x[SimData$train==TRUE],
           SimData$yobs[SimData$train==TRUE], sp, 1))-yobs[SimData$train==FALSE])^2))
ggplot(data.frame(span = spans, SSE = SSEfull), aes(x = span, y = SSE)) + geom_line()


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
spans[which.min(SSEfull)]


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
p + geom_line(data = SmoothData[SmoothData$span==spans[which.min(SSEfull)], ],
              aes(x = grid, y = value), color = "red")


## A regresszió legtöbb alkalmazott statisztikai terület talán legfontosabb eszköze

## 
## **Regresszió**: változók közti kapcsolat (illetve annak becslése minta alapján)

## 
## ,,Kapcsolat" formalizálása: függvény a matematikai fogalmával, tehát keressük az

## \[

## Y=f\left(X_1,X_2,\ldots,X_p\right)+\varepsilon=f\left(\mathbf{X}\right)

## \]

## függvényt

## 
## ($Y$ eredményváltozó, $X_i$-k a magyarázó változók)


## **Paraméteres regresszió**: ha *a priori* feltételezzük, hogy az $f$ függvény valamilyen -- paraméterek erejéig meghatározott -- függvényformájú (az ,,alakja" ismert), és így a feladat e paraméterek becslésére redukálódik

## 
## Tipikus példa a **lineáris regresszió**: $f\left(\mathbf{X}\right)=\beta_0+\beta_1 X_1 + \beta_2 X_2 + \ldots + \beta_p X_p=\mathbf{X}^T\pmb{\beta}$, így $Y=\mathbf{X}^T\pmb{\beta}+\varepsilon$

## 
## Ha rendelkezésre állnak az $\left\{y_i,\mathbf{x}_i\right\}_{i=1}^n$ megfigyeléseink a háttéreloszlásra, akkor e mintából megbecsülhetjük a paramétereket például **hagyományos legkisebb négyzetek** (OLS) módszerével:

## \[

## \widehat{\pmb{\beta}}=\argmin_{\mathbf{b}} \sum_{i=1}^n \left[Y_i-\mathbf{X}_i^T\mathbf{b}\right]^2=\left\| \mathbf{Y} - \mathbf{X}\mathbf{b} \right\|^2

## \]

## 
## Itt tehát $\mathbf{X}$ az a mátrix, amiben a magyarázó változók elé egy csupa 1 oszlopot szúrtunk, a neve **modellmátrix** vagy design mátrix


## De cserében mindig ott lebeg felettünk a kérdés, hogy a függvényformára *jó feltételezést* tettünk-e (hiszen ez nem az adatokból következik, ezt ,,ráerőszakoljuk" az adatokra)

## 
## (Persze ezért van a modelldiagnosztika)

## 
## A nem-paraméteres regresszió *flexibilis*, olyan értelemben, hogy minden a priori megkötés nélkül követi azt, ami az adatokból következik (a valóság ritkán lineáris?)

## 
## Cserében nehezebb becsülni, és nem kapunk analitikus -- jó esetben valamire hasznosítható -- regressziós függvényt, nem lehet értelmesen interpolálni és extrapolálni (,,fordul a kocka" a paraméteres esethez képest)


## Maradva a paraméteres keretben, arra azért mód van, hogy a függvényformát kibővítsük (és így flexibilisebbé tegyük)

## 
## Ezzel a különféle **nemlineáris regressziókhoz** jutunk el

## 
## E nemlinearitásoknak két alaptípusa van

## 
## -   Változójában nemlineáris modell (pl. $\beta_0 + \beta_1 x + \beta_2 x^2$): csak a szó ,,matematikai értelmében" nemlineáris, ugyanúgy becsülhető OLS-sel

## 
## -   Paraméterében nemlineáris modell (pl. $\beta_0x_1^{\beta_1}x_2^{\beta_2}$): felrúgja a lineáris struktúrát, így érdemileg más, csak linearizálás után, vagy NLS-sel becsülhető

## 
## Mi most az első esettel fogunk foglalkozni

## 
## Az itt látott ,,polinomiális regresszió" valóban nagyon gyakori módszer a flexibilitás növelésére


## Tekintsünk most egy másik példát, egy zajos másodfokú függvényt, kevesebb pontból:


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
n <- 20
x <- runif(n, 0, 10)
xgrid <- seq(0, 10, length.out = 100)
ygrid <- xgrid^2
yobs <- x^2 + rnorm(n, 0, 5)
SimData <- data.frame(x, xgrid, ygrid, yobs)
p <- ggplot(SimData) + geom_point(aes(x = x, y = yobs)) +
  geom_line(aes(x = xgrid, y = ygrid), color = "orange", lwd = 1)
p


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
fit5 <- lm(yobs ~ poly(x, 5), data = SimData)
p + geom_line(data = data.frame(xgrid, pred = predict(fit5, data.frame(x = xgrid))),
              aes(x = xgrid, y = pred))


## Mondjuk, hogy nagyobb flexibilitásra vágyunk

## 
## -   Például figyelembe akarjuk venni, hogy ez nem tűnik teljesen lineárisnak, vagy meg akarjuk ragadni a finomabb tendenciákat is

## 
## Emeljük a polinom fokszámát (ez nyilván növeli a flexibilitást, hiszen a kisebb fokszám nyilván speciális eset lesz), például 10-re

## 
## Szokás azt mondani, hogy a rang 5 illetve 10 (a polinom fokszáma, a becsülendő paraméterek száma nyilván egyezik a modellmátrix rangjával, de ez a fogalom később, amikor nem is polinomunk van, akkor is használható)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
fit10 <- lm(yobs ~ poly(x, 10), data = SimData)
p + geom_line(data = data.frame(xgrid, pred = predict(fit10, data.frame(x = xgrid))),
              aes(x = xgrid, y = pred))


## Szokás azt mondani, hogy *túlilleszkedés*, ami persze igaz is, de itt többről van szó

## 
## A polinomok elsősorban *lokálisan* tudnak jól közelíteni (a Taylor-sorfejtéses érvelés miatt), de nekünk arra lenne szükségünk, hogy *globálisan* jól viselkedő függvényformát találjunk

## 
## Pedig a polinomokat amúgy szeretjük, többek között azért is, mert szép sima görbét írnak le (matematikai értelemben véve a simaságot: végtelenszer folytonosan deriválhatóak, $C^{\infty}$-beliek)

## 
## Mi lehet akkor a megoldás?


## Egy lehetséges megközelítés: ,,összerakjuk a globálisat több lokálisból"

## 
## Azaz szakaszokra bontjuk a teljes intervallumot, és mindegyiket *külön-külön* polinommal igyekszünk modellezni

## 
## Így próbáljuk kombinálni a két módszer előnyeit

## 
## Persze a szakaszosan definiált polinomok önmagában még nem jók: a szakaszhatárokon találkozniuk kell (e találkozópontok neve: **knot**, ,,csomópont", a számukat $q-2$-val jelöljük, a pozíciójukat $x_i^{\ast}$-vel)

## 
## Sőt, ha a simasági tulajdonságokat is át akarjuk vinni, akkor az érintkezési pontokban a deriváltaknak (magasabbrendűeknek is) is egyezniük kell

## 
## Ha $p$-edfokú polinomokat használunk, akkor az első $p-1$ derivált -- és persze a függvényérték -- egyezését kell kikötnünk a knot-okban (és esetleg még valamit a végpontokra)

## 
## Ez így már jó konstrukció lesz, a neve: **spline**


## (Azért köbös, mert harmadfokúak a polinomok, és azért természetes, mert azt kötöttük ki, hogy a végpontokban nulla legyen a második derivált)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
fitSpline <- lm(yobs ~ splines::ns(x, 10), data = SimData)
p + geom_line(data = data.frame(xgrid, pred = predict(fitSpline, data.frame(x = xgrid))),
              aes(x = xgrid, y = pred))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
p + geom_line(data = rbind(data.frame(type = "Ötödfokú polinom",
                                      pred = predict(fit5, data.frame(x = xgrid)), xgrid),
                           data.frame(type = "Tizedfokú polinom",
                                      pred = predict(fit10, data.frame(x = xgrid)), xgrid),
                           data.frame(type = "Spline",
                                      pred = predict(fitSpline, data.frame(x = xgrid)),
                                      xgrid)),
              aes(x = xgrid, y = pred, color = type)) + labs(color = "")


## Nem csak az a jó, hogy szépen illeszkedik (tulajdonképpen még annál is jobban, mint a tizedfokú polinom, még ott is, ahol az jól illeszkedik amúgy)

## 
## ...hanem, hogy -- most már elárulhatom -- *ez is ugyanúgy 10 rangú* mint a tizedfokú polinom!

## 
## Mégis: nyoma nincs túlilleszkedésnek


## Amiről nem beszéltünk eddig: ez mind szép, de hogyan tudunk ténylegesen is megbecsülni egy ilyen spline-regressziót?

## 
## Ehhez visszalépünk pár lépést, és bevezetünk egy első kicsit absztraktnak tűnő, de később rendkívül jó szolgálatot tevő megközelítést

## 
## Bár a célunk a spline-regresszió becslésének a megoldása, de a dolog -- értelemszerűen -- alkalmazható polinomiális regresszióra is (legfeljebb nincs sok értelme, mert az hagyományos módszerekkel is jól kézbentartható), úgyhogy először azon fogjuk illusztrálni


## A másodfokú polinomok -- mint függvények -- összessége **függvényteret** alkot

## 
## Ez egy olyan *vektortér*, aminek az elemei a függvények, a skalárok a valós számok, a két művelet pedig

## 
## -   Skalárral szorzás: $\left(cf\right)\left(x\right)=cf\left(x\right)$

## 
## -   Vektorok (azaz függvények) összeadása: $\left(f+g\right)\left(x\right)=f\left(x\right)+g\left(x\right)$, tehát pontonkénti összeadás

## 
## Belátható, hogy ez teljesíti a vektortéraxiómákat, mert zárt a két műveletre (másodfokú polinomok összege másodfokú polinom és másodfokú polinom konstansszorosa másodfokú polinom), és a többi követelményt is teljesíti


## Szuper, de mindez mire jó?

## 
## Ha vektortér, akkor létezik **bázisa**, azaz olyan vektorok halmaza, melyekből lineáris kombinációval minden vektor -- egyértelműen -- előállítható (bázis: lineárisan független generátorrendszer)

## 
## A bázis nem feltétlenül egyértelmű, de az elemszáma igen, ez a vektortér **dimenziója**

## 
## Például a másodfokú polinomok jó bázisa $\left\{1,x,x^2\right\}$, nyilvánvaló, hogy ebből tényleg minden $ax^2+bx+c$ másodfokú polinom előállítható lineáris kombinációval (triviálisan, a súlyok $c$, $b$ és $a$)

## 
## Függvényterek esetében a bázis elemeit **bázisfüggvényeknek** is szokás nevezni, az $\left\{1,x,x^2\right\}$ tehát a másodfokú polinomok bázisfüggvényei


## Mivel mutattunk egy konkrét bázist, így a dimenzió nyilván 3, de a későbbiek szempontjából jól jön egy másik módszer is

## 
## Azzal, hogy az $ax^2+bx+c$ polinomot megfeleltettük az $\left(a,b,c\right)$ valós számhármasnak, a polinomok tere és a valós számhármasok tere (az $\mathbb{R}^3$) között létesítettünk egy izomorfizmust (a leképezés művelettartó és kölcsönösen egyértelmű)

## 
## Emiatt a polinomok terének ugyanaz a dimenziója, mint az $\mathbb{R}^3$-nak, ami viszont természetesen 3

## 
## Ez a módszer általában is használható: a dimenzió a felíráshoz szükséges paraméterek száma (feltéve, hogy ezek valós számok, valamint mindegyikhez tartozik egy polinom és viszont)


## Mindez a spline-okra is igaz!

## 
## Érthető: minden pontban két polinomot adunk össze, vagy polinomot szorzunk skalárral, az eredmény polinom (már láttuk) -- így tud spline adott pontja lenni!

## 
## Azaz: spline-okat is elő tudunk állítani bázisfüggvények lineáris kombinációjaként!


## Mielőtt megkeressük a spline-ok terének egy bázisát (azaz a konkrét bázisfüggvényeket), tisztázni kellene, hogy hány bázisfüggvényt keresünk egyáltalán, azaz hány dimenziós a spline-ok függvénytere

## 
## Naiv ötlet (köbös spline-okat használva példaként): van $q-1$ szakasz ($q-2$ knot, ami meghatároz $q-3$ szakaszt meg a két vége; úgy is felfogható, hogy a két végével együtt $q$ knot van, ami meghatároz $q-1$ szakaszt) és mindegyiken egy harmadfokú polinom (aminek 4 paramétere van), akkor az $4q-4$ paraméter

## 
## Igen ám, de vannak megkötések: a knotokban a függvényérték és az első két derivált egyezik

## 
## Minden megkötés minden pontban 1 egyenlet, az 1-gyel csökkenti a paraméterek számát: van $q-2$ knot és 3 megkötés, az $3q-6$ csökkentés, marad $q+2$ paraméter

## 
## De mivel természetes, így a végpontokban is van 1-1 megkötés: marad $q$ paraméter, azaz $q$ dimenziós a természetes köbös spline-ok tere (ezért neveztük a knot-ok számát $q-2$-nek!)


## Természetesen itt is igaz, hogy adott, rögzített spline-osztályra (pl. természetes köbös) is végtelen sok bázis van

## 
## Köztük célszerűség alapján választhatunk

## 
## A részletek nélkül két példa:

## 
## -   $b_1\left(x\right)=1, b_2\left(x\right)=x, b_i\left(x\right)=\left|x-x_{i-2}^{\ast}\right|^3 (i=3,4,\ldots,q)$

## 
## -   $b_1\left(x\right)=1, b_2\left(x\right)=x, b_i\left(x\right)=R\left(x,x_{i-2}^{\ast}\right) (i=3,4,\ldots,q)$, ahol $R$ egy nevezetes -- elég hosszú, bár nem túl bonyolult -- függvény (hamar látni fogjuk, hogy ez miért előnyös), annyi fontos, hogy $x$ a $\left[0,1\right]$ intervallumban essen (egyszerű átskálázssal mindig elérhető)

## 
## Most már csak a regresszió kivitelezését kell kitalálnunk


## A bázisfüggvények használatának két hatalmas előnye van:

## 
## -   A probléma visszavezethető velük a sima lineáris regresszióra

## 
## -   Sőt, ehhez a modellmátrix is könnyen előállítható


## Legyen $b_1\left(x\right)=1$, $b_2\left(x\right)=x$ és $b_3\left(x\right)=x^2$ a bázisunk

## 
## Az eredeti regresszió:

## \[

## y_i = \beta_1 + \beta_2 x_i + \beta_3 x_i^2 + \varepsilon_i

## \]

## 
## Átírva bázisokra (lényegében transzformált magyarázó változók):

## \[

## y_i = \beta_1 b_1\left(x_i\right) + \beta_2 b_2\left(x_i\right) + \beta_3 b_3\left(x_i\right) + \varepsilon_i

## \]

## 
## Ez már tiszta lineáris regresszió


## Ez úgy tűnik, hogy csak egy nagyon nyakatekert felírás egy amúgy egyszerű problémára

## 
## Valójában viszont egy elképesztően erőteljes dolgot nyertünk: *minden* olyan függvény, legyen bármilyen komplikált is, ami felírható bázisfüggvényekkel (azaz az osztálya függvényosztályt alkot), az berakható egy *kutyaközönséges* regresszióba (azaz lehet ő a regrssziós függvény) a fenti transzformációval, tehát

## \[

##     \sum_{i=1}^q \beta_i b_i\left(x\right)

## \]

## alakban

## 
## (Azaz minden függvény, ami egy függvénytér eleme)


## Még egyszer: *minden* függvény, ami felírható bázisfüggvényekkel

## 
## Azaz: *minden*

## 
## ...és az összesnek *pontosan ugyanúgy* az lesz az alakja, hogy

## \[

## \sum_{i=1}^q \beta_i b_i\left(x\right),

## \]

## egyedül a bázisfüggvényt kell az adott esetnek megfelelően megválasztani

## 
## Tehát a spline is mehet ugyanígy (csak megfelelő $b_i$-kkel)!

## 
## És ha ez az alak megvan, akkor onnantól természetesen *sima lineáris regresszióval* elintézhető


## Ráadásul az $\mathbf{X}$ modellmátrix (design mátrix) előállítása is nagyon könnyű lesz: az $i$-edik sora

## \[

## \left[ b_1\left(x_i\right), b_2\left(x_i\right), \ldots, b_q\left(x_i\right) \right]

## \]

## 
## Így maga a mátrix az $\mathbf{x}$ és az $\left[1,2,\ldots,q\right]$ vektor *külső szorzata* (tenzorszorzata), ha a művelet alatt az oszlopban szereplő érték által meghatározott bázisfüggvény sorbeli elemre történő alkalmazását értjük, tehát $i\otimes j:=b_j\left(x_i\right)$, és így

## \[

## \begin{aligned}

## & \begin{pmatrix} \quad 1 & \qquad \enspace 2 & \quad \; \cdots & \quad q\quad \; \end{pmatrix} \\

## \begin{pmatrix}x_1\\x_2\\ \vdots \\ x_n\end{pmatrix} & \begin{bmatrix}b_1\left(x_1\right) & b_2\left(x_1\right) & \cdots & b_q\left(x_1\right) \\ b_1\left(x_2\right) & b_2\left(x_2\right) & \cdots & b_q\left(x_2\right) \\ \vdots & \vdots & \ddots & \vdots \\ b_1\left(x_n\right) & b_2\left(x_n\right) & \cdots & b_q\left(x_n\right) \end{bmatrix}

## \end{aligned}

## \]

## 
## Így, a teljes modellmátrix egy lépésben megkapható...

## 
## ... majd közvetlenül rakható is bele a sima lineáris regresszióba (ld. 1. előny):

## \[

## \widehat{\boldsymbol{\beta}}=\left(\mathbf{X}^T\mathbf{X}\right)^{-1}\mathbf{X}^T\mathbf{y}

## \]


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
n <- 30
x <- runif(n, 0, 1)
xgrid <- seq(0, 1, length.out = 100)
ygrid <- 100*xgrid^2
yobs <- 100*x^2 + rnorm(n, 0, 5)
p <- ggplot(data.frame(x, yobs)) + geom_point(aes(x = x, y = yobs)) +
    geom_line(data = data.frame(xgrid, ygrid), aes(x = xgrid, y = ygrid),
              color = "orange", lwd = 1)
p


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
xk <- 1:4/5
q <- length(xk) + 2


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
rk <- function( x, z ) {
    ((z-0.5)^2-1/12)*((x-0.5)^2-1/12)/4-((abs(x-z)-0.5)^4-(abs(x-z)-0.5)^2/2+7/240)/24
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
X <- matrix(1, n, q) 


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
X[, 2] <- x


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
X[, 3:q] <- outer(x, xk, FUN = rk)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
spl.X <- function(x, xk) {
    q <- length(xk) + 2
    n <- length(x)
    X <- matrix(1, n, q)
    X[, 2] <- x
    X[, 3:q] <- outer(x, xk, FUN = rk)
    X
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
fit <- lm(yobs ~ X - 1 )


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
Xp <- spl.X(xgrid, xk)
yp <- Xp%*%coef(fit)
p + geom_line(data = data.frame(xgrid, yp), aes(x = xgrid, y = yp))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
predspline <- function(x, y, q) {
    xk <- (1:(q-2))/(q-1)
    X <- spl.X(x, xk)
    fit <- lm(y ~ X - 1)
    xp <- 0:100/100
    Xp <- spl.X(xp, xk)
    yp <- Xp%*%coef(fit)
    list(fit = fit, xp = xp, yp = yp)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
p + geom_line(data = with(predspline(x, yobs, 6), data.frame(xp, yp)),
              aes(x = xp, y = yp))
p + geom_line(data = with(predspline(x, yobs, 11), data.frame(xp, yp)),
              aes(x = xp, y = yp))
p + geom_line(data = with(predspline(x, yobs, 3), data.frame(xp, yp)),
              aes(x = xp, y = yp))


## A $q$ dimenzió tehát az illeszkedés szabadságát határozza meg

## 
## Valahogy ezt is meg kellene határozni

## 
## Jön a fő kérdéskör: a túlilleszkedés elleni védekezés

## 
## Milyen legyen a ,,simítás foka"?


## Tehát $q$-t kellene valahogy jól belőni

## 
## Egyszerű *modellszelekció*?

## 
## -   Vagy nem beágyazott modellek szelekciója, vagy nem ekvidisztáns knot-ok, egyik sem túl szerencsés

## 
## Alternatív ötlet: $q$ legyen inkább rögzített (elég nagy értéken, kicsit a várható fölé lőve), de a függvényformát nem engedjük teljesen szabadon alakulni

## 
## Hogyan? Büntetjük a túl ,,zizegős" függvényt!

## 
## Ez épp a **penalizált regresszió** alapötlete

## 
## És ami rendkívül fontos: így már jellemzően sem $q$ pontos megválasztása, sem a knot-ok pontos helye nem bír nagy jelentőséggel (választhatjuk például egyenletesen)!


## Klasszikus megoldás: a második derivált jelzi adott pontban a ,,zizegősséget", ezt kiintegrálva kapunk egy összesített mértéket az egész függvényre

## 
## Valamilyen súllyal ezt vegyük figyelembe:

## \[

##     \left\|\mathbf{y}-\mathbf{X}\boldsymbol{\beta}\right\|^2+\lambda\int_0^1 \left[f''\left(x\right)\right]^2 \dif x

## \]

## 
## A $\lambda$ a *simítási paraméter*, ez határozza meg a trade-off-ot a jó illeszkedés és a simaság között

## 
## -   $\lambda=0$: penalizálatlan becslés, $\lambda\rightarrow\infty$: egyenes regressziós függvény


## A regressziós függvény alakja: $f\left(x\right)=\sum_{i=1}^q \beta_i b_i\left(x\right)$

## 
## Kétszer deriválva: $f''\left(x\right)=\sum_{i=1}^q \beta_i b_i''\left(x\right)$

## 
## Négyzetre emelve: $\left[f''\left(x\right)\right]^2=\sum_{i=1}^q \sum_{j=1}^q \beta_i b_i''\left(x\right) b_j''\left(x\right) \beta_j$

## 
## Kiintegrálva: $\int_0^1 \left[f''\left(x\right)\right]^2 \dif x=\sum_{i=1}^q \sum_{j=1}^q \beta_i \left(\int_0^1 b_i''\left(x\right) b_j''\left(x\right) \dif x\right) \beta_j$

## 
## De hát ez épp egy *kvadratikus alak*! ($\sum_{i=1}^q \sum_{j=1}^q x_i a_{ij} x_j= \mathbf{x}^T \mathbf{A} \mathbf{x}$)

## 
## Legyen $S_{ij}=\int_0^1 b_i''\left(x\right) b_j''\left(x\right) \dif x$ és $\mathbf{S}$ az ezekből alkotott mátrix, akkor tehát a simítási büntetőtag:

## \[

##     \lambda \boldsymbol{\beta}^T\mathbf{S}\boldsymbol{\beta}

## \]

## 
## Az előbb definiált $R$-rel $\mathbf{S}$ alakja nagyon egyszerű lesz: $S_{i+2,j+2}=R\left(x_i^{\ast},x_j^{\ast}\right)$, az első két oszlop és sor pedig csupa nulla


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
spl.S <- function(xk) {
    q <- length(xk) + 2
    S <- matrix(0, q, q)
    S[3:q, 3:q] <- outer(xk, xk, FUN = rk)
    S
}


## Kényelmes lenne, ha $\left\|\mathbf{y}-\mathbf{X}\boldsymbol{\beta}\right\|^2+\lambda \boldsymbol{\beta}^T\mathbf{S}\boldsymbol{\beta}$ helyett írhatnánk egyetlen normát célfüggvényként

## 
## Ez nem nehéz, ha a második tagot át tudjuk normává alakítani, hiszen (innentől némi blokkmátrix műveletekre szükség lesz)

## \[

## \left\|\mathbf{a}\right\|^2+\left\|\mathbf{b}\right\|^2=\left\|\begin{pmatrix}\mathbf{a} \\ \mathbf{b} \end{pmatrix}\right\|^2

## \]

## 
## Legyen $\mathbf{B}$ olyan, hogy $\mathbf{B}^T\mathbf{B}=\mathbf{S}$ (pl. spektrális dekompozícióval, vagy Cholesky-dekompozícióval megtalálható a mátrix ilyen ,,négyzetgyöke"), ekkor

## \[

## \lambda \boldsymbol{\beta}^T\mathbf{S}\boldsymbol{\beta} = \lambda \boldsymbol{\beta}^T\mathbf{B}^T\mathbf{B}\boldsymbol{\beta}=\lambda \left( \mathbf{B} \boldsymbol{\beta}\right)^T\mathbf{B}\boldsymbol{\beta} =\left( \sqrt{\lambda} \mathbf{B} \boldsymbol{\beta}\right)^T\left(\sqrt{\lambda}\mathbf{B}\boldsymbol{\beta}\right)

## \]

## 
## Ezzel meg is vagyunk, hiszen a norma egyszerűen $\left\|\mathbf{a}\right\|^2=\mathbf{a}^T\mathbf{a}$, így

## \[

## \lambda \boldsymbol{\beta}^T\mathbf{S}\boldsymbol{\beta} = \left\|\sqrt{\lambda} \mathbf{B} \boldsymbol{\beta}\right\|^2

## \]

## ahonnan

## \[

## \left\|\mathbf{y}-\mathbf{X}\boldsymbol{\beta}\right\|^2+\lambda \boldsymbol{\beta}^T\mathbf{S}\boldsymbol{\beta}=\left\|\mathbf{y}-\mathbf{X}\boldsymbol{\beta}\right\|^2+\left\|\sqrt{\lambda} \mathbf{B} \boldsymbol{\beta}\right\|^2\]

## és így, az előzőek szerint

## \[

## \left\|\mathbf{y}-\mathbf{X}\boldsymbol{\beta}\right\|^2+\left\|\sqrt{\lambda} \mathbf{B} \boldsymbol{\beta}\right\|^2=\left\|\begin{pmatrix}\mathbf{y}-\mathbf{X}\boldsymbol{\beta} \\ \sqrt{\lambda} \mathbf{B} \boldsymbol{\beta} \end{pmatrix}\right\|^2

## \]

## 
## Jó lenne $\boldsymbol{\beta}$-t kiemelni; ez nem is túl nehéz, hiszen $\mathbf{a}$ és $-\mathbf{a}$ normája ugyanaz:

## \[

## \left\|\begin{pmatrix}\mathbf{y}-\mathbf{X}\boldsymbol{\beta} \\ \sqrt{\lambda} \mathbf{B} \boldsymbol{\beta} \end{pmatrix}\right\|^2 = \left\| \begin{pmatrix} \mathbf{y} \\ \mathbf{0} \end{pmatrix} - \begin{pmatrix}\mathbf{X} \\ \sqrt{\lambda} \mathbf{B} \end{pmatrix}\boldsymbol{\beta}\right\|^2

## \]


## Innentől a regresszió játszi könnyedséggel (értsd: a szokványos, nem is penalizált eszköztárral) megoldható, csak $\mathbf{X}$ szerepét $\begin{pmatrix}\mathbf{X} \\ \sqrt{\lambda} \mathbf{B} \end{pmatrix}$, $\mathbf{y}$ szerepét $\begin{pmatrix} \mathbf{y} \\ \mathbf{0} \end{pmatrix}$ játssza

## 
## Így az ,,$\mathbf{X}^T\mathbf{X}$" épp $\mathbf{X}^T\mathbf{X}+\lambda \mathbf{B}^T\mathbf{B}=\mathbf{X}^T\mathbf{X}+\lambda\mathbf{S}$ lesz

## 
## Az ,,$\mathbf{X}^T\mathbf{y}$" pedig $\mathbf{X}^T\mathbf{y}$ (a kiegészített eredményváltozóban lévő nullák épp a magyarázó változók kiegészítését ütik ki)

## 
## Így az OLS megoldás:

## \[

## \widehat{\boldsymbol{\beta}}=\left(\mathbf{X}^T\mathbf{X}+\lambda\mathbf{S}\right)^{-1}\mathbf{X}^T\mathbf{y}

## \]

## 
## (Persze a gyakorlatban ennek közvetlen számítása helyett célszerűbb az augmentált eredmény- és magyarázóváltozókat berakni egy hatékonyabb lineáris regressziót megoldó módszerbe)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
mat.sqrt <- function(S) {
    d <- eigen(S, symmetric = TRUE)
    d$vectors%*%diag(d$values^0.5)%*%t(d$vectors)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
predsplinepen <- function(x, y, q, lambda) {
    xk <- (1:(q-2))/(q-1)
    Xa <- rbind(spl.X(x, xk), sqrt(lambda) * mat.sqrt(spl.S(xk)))
    ya <- c(y, rep(0, q))
    fit <- lm(ya ~ Xa - 1)
    xp <- 0:100/100
    Xp <- spl.X(xp, xk)
    yp <- Xp%*%coef(fit)
    list(fit = fit, xp = xp, yp = yp)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
p + geom_line(data = with(predsplinepen(x, yobs, 20, 1), data.frame(xp, yp)),
              aes(x = xp, y = yp))
p + geom_line(data = with(predsplinepen(x, yobs, 20, 0.001), data.frame(xp, yp)),
              aes(x = xp, y = yp))
p + geom_line(data = with(predsplinepen(x, yobs, 20, 0.000001), data.frame(xp, yp)),
              aes(x = xp, y = yp))


## Kérdés még a $\lambda$ értéke

## 
## Sima OLS-jellegű eljárással, tehát a reziduális négyzetösszeg minimalizálást tűzve ki célul nyilván nem határozható meg (hiszen az mindig 0-t adna)

## 
## Épp az a lényeg, hogy a túlilleszkedésre is tekintettel legyünk

## 
## Ötlet: keresztvalidáció


## Mindig egy pontot hagyunk ki, és így számolunk hibát: OCV

## 
## (Szokták egy-kihagyásos keresztvalidációnak, LOOCV-nek is nevezni)

## 
## Tehát:

## \[

##     E_{OCV}=\frac{1}{n}\sum_{i=1}^n \left( \widehat{f}_i^{\left[-i\right]} - y_i\right)^2

## \]

## 
## Szerencsére nem kell ténylegesen $n$-szer lefuttatni a regressziót mert belátható, hogy

## \[

##     E_{OCV}=\frac{1}{n}\sum_{i=1}^n \left(y_i - \widehat{f}_i\right)^2/\left(1-A_{ii}\right)^2,

## \]

## ahol $\mathbf{A}$ az influence mátrix


## Ha az $A_{ii}$-ket az átlagukkal helyettesítjük, akkor az általánosított keresztvalidációhoz jutunk (GCV)

## 
## Tehát:

## \[

##     E_{GCV}=\frac{1}{n}\sum_{i=1}^n \left(y_i - \widehat{f}_i\right)^2/\left[\mathrm{tr}\left(\mathbf{I}-\mathbf{A}\right)\right]^2

## \]


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
predV <- 10^(seq(-8, 3, length.out = 100))
V <- sapply(predV, function(lambda) {
    fit <- predsplinepen(x, yobs, 20, lambda)$fit
    trA <- sum(influence(fit)$hat[1:n])
    rss <- sum((yobs - fitted(fit)[1:n])^2)
    n*rss/(n - trA)^2
} )
ggplot(data.frame(predV, V), aes(x = predV, y = V)) + geom_line() + scale_x_log10()


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
predV[which.min(V)]


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
p + geom_line(data = with(predsplinepen(x, yobs, 20, predV[which.min(V)]),
                          data.frame(xp, yp)), aes(x = xp, y = yp))


## Eddig egy magyarázó változó esetével foglalkoztunk

