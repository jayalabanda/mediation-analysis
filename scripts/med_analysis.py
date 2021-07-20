import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

from statsmodels.api.families import Binomial


def ind(x, value):
    res = [value if x[i] == value else 0 for i in range(len(x))]
    return res


def medtmle(data, covars, A, Z, M, outcome, amodel, zmodel, mmodel, ymodel, qmodel):
    obsdat = pd.copy(data)
    obsdat.rename(columns={A: "a", Z: "z", M: "m", outcome: "y"},
                  inplace=True)

    for i in range(1, len(covars) + 1):
        obsdat.rename(columns={covars[i]: "w" + str(i)},
                      inplace=True)

    tmpdat = obsdat.copy()
    dfm0 = tmpdat.copy()
    dfm1 = dfm0.copy()
    dfa0z0 = dfm1.copy()
    dfa0z1 = dfa0z0.copy()
    dfa1z0 = dfa0z1.copy()
    dfa1z1 = dfa1z0.copy()
    dfa0 = dfa1z1.copy()
    dfa1 = dfa0.copy()

    dfa1.a = 1
    dfa1z1.a = 1
    dfa1z1.z = 1
    dfa1z0.a = 1
    dfa0z1.z = 1
    dfm1.m = 1

    dfa0.a = 0
    dfa1z0.z = 0
    dfa0z1.a = 0
    dfa0z0.a = 0
    dfa0z0.z = 0
    dfm0.m = 0

    z_fit = smf.glm(zmodel, family=Binomial(), data=obsdat).fit()
    z_A0 = z_fit.predict(dfa0)
    z_A1 = z_fit.predict(dfa1)

    m_fit = smf.glm(mmodel, family=Binomial(), data=obsdat).fit()
    m_Z1A0 = m_fit.predict(dfa0z1)
    m_Z0A0 = m_fit.predict(dfa0z0)
    m_Z1A1 = m_fit.predict(dfa1z1)
    m_Z0A1 = m_fit.predict(dfa1z0)

    gm_A0 = m_Z1A0 * z_A0 + m_Z0A0 * (1 - z_A0)
    gm_A1 = m_Z1A1 * z_A1 + m_Z0A1 * (1 - z_A1)

    a_fit = smf.glm(amodel, family=Binomial(), data=tmpdat).fit()
    pred = a_fit.predict(a_fit)
    ps_A1 = ind(tmpdat.a, 1) / pred
    ps_A0 = ind(tmpdat.a, 1) / (1 - pred)

    m_fit = smf.glm(mmodel, family=Binomial(), data=tmpdat).fit()
    mazw = m_fit.predict()
    psm = ind(tmpdat.m, 1) * mazw + ind(tmpdat.m, 0) * (1 - mazw)

    y_fit = smf.glm(ymodel, family=Binomial(), data=tmpdat).fit()
    tmpdat['qyinit'] = tmpdat.append(
        [y_fit.predict(tmpdat), y_fit.predict(dfm0), y_fit.predict(dfm1)]
    )

    tmpdat['hA1gmA0'] = (
        (ind(tmpdat.m, 1) * gm_A0 + ind(tmpdat.m, 0) * (1 - gm_A0)) / psm) * ps_A1

    epsilon_A1gmA0 = smf.glm("y ~ 1", freq_weights=tmpdat.hA1gmA0, offset=)
