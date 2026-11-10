# Diagnostico del signo negativo de `B_T` en los modelos `droplet` y Diego

## Resumen

Se analizo por que los modelos `droplet` y `evaporacion_gota_diego.py` producen valores negativos de `B_T`, mientras que `evaporacion_gota_new.py` (`legacy_original`) produce valores positivos para el mismo caso.

La conclusion principal es:

- La ecuacion de `B_T` es esencialmente la misma en los tres modelos.
- La diferencia no nace en la formula de `B_T`.
- La diferencia nace en la temperatura superficial `Ts` que resuelve cada modelo.
- `legacy_original` cae en un valor de respaldo (`fallback`) para `Ts`.
- `droplet` y Diego si encuentran la raiz caliente del balance superficial.
- Esa raiz caliente hace que `Ts > T_bar`, lo que vuelve negativo a `q_l` y puede volver negativo a `B_T`.

## Formula de `B_T`

En los tres modelos se usa la misma estructura:

```text
B_T_den = DHev + q_l / m_dot
B_T     = cpg * (T_inf - Ts) / B_T_den
```

Referencias:

- [evaporacion_gota_new.py](/d:/programacion/evaporacion/python/src/evaporacion/evaporacion_gota_new.py#L257)
- [algebraic.py](/d:/programacion/evaporacion/python/src/evaporacion/droplet/algebraic.py#L157)
- [evaporacion_gota_diego.py](/d:/programacion/evaporacion/python/src/evaporacion/evaporacion_gota_diego.py#L325)

Por lo tanto, el signo de `B_T` depende de:

- el numerador `cpg * (T_inf - Ts)`
- el denominador `DHev + q_l / m_dot`

## Mecanismo que cambia el signo

El termino clave es el flujo de calor interno:

```text
q_l = -4 * pi * b_q * k_l * (Ts - T_bar) * rS
```

Si:

- `Ts < T_bar`, entonces `q_l > 0`
- `Ts > T_bar`, entonces `q_l < 0`

Como `m_dot > 0` durante evaporacion, cuando `q_l` es suficientemente negativo, el termino:

```text
q_l / m_dot
```

puede dominar a `DHev` y hacer negativo el denominador:

```text
DHev + q_l / m_dot < 0
```

En ese caso `B_T` se vuelve negativo.

## Diferencia real entre los modelos

### `legacy_original`

En `legacy_original`, la temperatura superficial se resuelve con un Brent sobre:

- intervalo inicial `[0.01, T_inf]`
- luego `[-10, T_inf + 5]` si no hay cambio de signo
- si falla, usa `Ts = 0.95 * T_bar`

Referencia:

- [evaporacion_gota_new.py](/d:/programacion/evaporacion/python/src/evaporacion/evaporacion_gota_new.py#L296)

El problema es que al evaluar cerca de `Ts = T_inf`, el modelo calcula:

- `B_T = 0`
- luego `ln(1 + B_T) / B_T`

sin regularizacion numerica.

Referencias:

- [evaporacion_gota_new.py](/d:/programacion/evaporacion/python/src/evaporacion/evaporacion_gota_new.py#L259)
- [evaporacion_gota_new.py](/d:/programacion/evaporacion/python/src/evaporacion/evaporacion_gota_new.py#L266)

Eso hace que el balance superficial en `Ts = T_inf` devuelva `NaN`, por lo que el solver termina cayendo en el `fallback`:

```text
Ts = 0.95 * T_bar
```

Para el caso diagnosticado, eso da:

```text
Ts = 20.9 C
```

Como `Ts < T_bar`, entonces:

- `q_l > 0`
- `DHev + q_l/m_dot > 0`
- `B_T > 0`

### `droplet`

En `droplet`, el balance superficial usa la misma fisica, pero introduce una regularizacion numerica de `B_T` antes de usar `ln(1 + B_T) / B_T`:

```text
_BT = B_T si |B_T| > 1e-12, o +/- 1e-12 si no
```

Referencia:

- [algebraic.py](/d:/programacion/evaporacion/python/src/evaporacion/droplet/algebraic.py#L160)

Eso evita el `NaN` en la evaluacion del balance y permite que el solver encuentre la raiz fisica caliente:

```text
Ts ~= 22.381 C
```

Como `Ts > T_bar`, entonces:

- `q_l < 0`
- `q_l / m_dot` vuelve negativo al denominador
- `B_T < 0`

### Diego

El modelo Diego mantiene la misma estructura de `B_T`, pero ademas permite resolver `Ts` con un limite superior mas alto:

```text
Ts_hi = T_inf + 10
```

Referencia:

- [evaporacion_gota_diego.py](/d:/programacion/evaporacion/python/src/evaporacion/evaporacion_gota_diego.py#L365)

Eso facilita encontrar una raiz con `Ts > T_bar`, evitando el comportamiento del `fallback` de `legacy_original`.

El resultado es el mismo patron observado en `droplet`:

- `Ts > T_bar`
- `q_l < 0`
- `B_T < 0`

## Evidencia del caso diagnosticado

Caso usado:

- `d0 = 2.8 mm`
- `T0 = 22 C`
- `T_inf = 35 C`
- `HR = 0.3`
- `U = 1 m/s`
- `perfil = parabolico`

Valores comparados en el estado inicial:

| Variable | `legacy_original` | `droplet` |
| --- | ---: | ---: |
| `Ts [C]` | `20.9` | `22.3811` |
| `q_l [W]` | `+0.1595` | `-0.05527` |
| `m_dot [kg/s]` | `9.5967e-09` | `1.24534e-08` |
| `DHev [J/kg]` | `2.45234e+06` | `2.44873e+06` |
| `B_T` | `+7.44e-04` | `-6.39e-03` |
| `B_M` | `5.10e-03` | `6.65e-03` |
| `Y_s` | `1.576e-02` | `1.727e-02` |
| `DHev + q_l/m_dot` | `+1.9076e+07` | `-1.9894e+06` |

Interpretacion:

- En `legacy_original`, el `fallback` fuerza una `Ts` mas baja.
- Esa `Ts` baja hace que `q_l` sea positivo.
- El denominador de `B_T` queda positivo.
- Por eso `B_T` sale positivo.

En cambio:

- En `droplet` y Diego, la raiz encontrada tiene `Ts` mayor que `T_bar`.
- Eso vuelve `q_l` negativo.
- El termino `q_l / m_dot` domina y vuelve negativo el denominador.
- Por eso `B_T` sale negativo.

## Confirmacion directa del problema numerico

Se verifico que el balance superficial en `legacy_original` da:

- valor finito en `Ts = 0.01`
- `NaN` en `Ts = T_inf = 35`
- valor finito en `Ts = 22.381...`

Eso demuestra que:

- si `legacy_original` pudiera atravesar sin problemas el punto singular en `Ts = T_inf`, tambien podria encontrar la raiz caliente
- el signo distinto de `B_T` no es una diferencia fisica primaria del modelo
- es principalmente una diferencia numerica en la resolucion de `Ts`

## Conclusion

La causa principal de los `B_T` negativos en `droplet` y Diego es que esos modelos si resuelven la rama caliente de la temperatura superficial, mientras que `legacy_original` cae en un `fallback` mas frio debido a una singularidad numerica cuando `B_T -> 0` en `Ts = T_inf`.

En una frase:

> `B_T` es negativo en `droplet` y Diego porque alli `Ts > T_bar`, lo que hace `q_l < 0` y vuelve negativo el denominador `DHev + q_l/m_dot`; `legacy_original` evita eso solo porque su solver de `Ts` cae en un valor de respaldo.

## Implicacion practica

Si se regulariza `legacy_original` del mismo modo que `droplet` al evaluar `B_T` cerca de cero, es esperable que:

- deje de caer en el `fallback`
- encuentre la misma raiz caliente
- y pase a producir `B_T` negativos o muy cercanos a los de `droplet`

