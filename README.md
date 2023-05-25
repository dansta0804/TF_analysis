# Transkripcijos faktorių taikinių spėjimo metodo kūrimas ir implementacija *Shiny* aplikacijoje

Repozitorijoje saugomi failai, kurie yra naudojami sukurtoje
R aplikacijoje. Ši aplikacija leidžia įvertinti įkeltų ChIP
sekoskaitos duomenų kokybę, atlikti biologines analizes bei
panaudoti sukurtą metodą transkripcijos faktorių taikinių
spėjimui pasirinktame organizme.

**Direktorijų apibūdinimai:**
- **Input/** - šioje direktorijoje saugomas JSON formato failas,
kuriame yra pateikta visa R bibliotekų, genomų anotacijų bei
chromosomų rinkinio ir jų ilgių informacija kiekvienam
organizmui, kurio ChIP sekoskaitos duomenys gali būti apdoroti
aplikacijoje.
- **Scripts/** - pagrindinė direktorija, kurioje saugomi visi
failai, realizuojantys sukurtą aplikaciją. Šioje direktorijoje
patalpintas aplikacijoje naudojamų funkcijų failas (*functions.R*),
aplikaciją sudaranačių elementų stilių apibūdinantis CSS formato
failas (*styles.css*) bei **App** direktorija. Šią direktoriją
sudaro trys R failai:
  - *ui.R* failas aprašytas aplikacijos elementų išdėstymas;
  - *server.R* faile realizuotas R programinis kodas,
  atliekantis ChIP sekoskaitos duomenų kokybės vertinimą,
  biologines analizes bei realizuojantis transkripcijos faktorių
  taikinių spėjimo metodą. Tai pat šiame faile esantis kodas
  generuoja iliustracijas, tarpinius programos failus bei
  lenteles;
  - *app.R* - pagrindinis failas, kuriame inicijuojamas
  *ui.R* ir *server.R* programų kvietimas. Šis failas reikalingas
  *Shiny* aplikacijos paleidimui.

**Aplikacijos paleidimas**
1. Norint paleisti sukurtą *Shiny* aplikaciją reikia atsisiųsti
virtualios mašinos **Aplikacija.OVA** failą, kuris gali būti atsisiųtas
iš čia:
[https://www.google.com](https://www.google.com).

2. Atsisiųstas failas turi būti įkeltas į Oracle VirtualBox:
**File** > **Import appliance** > **File** dalyje reikia
įkelti atsisiųstą **Aplikacija.OVA** failą.

3. Pasirodžius darbalaukiui reikia pasirinkti aplanką "Aplikacija".

4. Šiame aplanke atverti terminalą ir įrašyti komandą:
`Rscript ./Scripts/App/app.R`

5. Pasirodžius *http* adresui jį įvesti į naudojamą naršyklę.
