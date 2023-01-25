import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
from tqdm import tqdm
import seaborn as sns
sns.set()

def La(x, a):
    return np.exp(a*x) - a*x - 1

x = np.linspace(-3,3,100)
for a in [0.1, 0.5, 1.2]:
    plt.plot(x, La(x, a), label="a="+str(a))
    plt.legend()
plt.title("Cout La pour differents a")
plt.show()

def rejection_sampling(likelihood, prior, instr, c, n_iter=100000, n_sample=100000):
    """
    :param likelihood: fct theta -> log( f(x|theta) )
    :param prior: fct theta -> pi(theta)
    :param instr: fct theta -> rho(theta)
    :param c:  log (sup theta -> f*pi/rho)
    """

    accepted_values = []
    acceptance_rate = []
    for i in tqdm(range(n_iter)):
        theta = instr.rvs()
        log_prob = likelihood(theta) + prior.logpdf(theta) - instr.logpdf(theta) - c
        #èprint(theta, log_prob)
        if np.log(np.random.uniform()) < log_prob:
            # accepted
            accepted_values.append(theta)
            acceptance_rate.append(1)
        else:
            acceptance_rate.append(0)
        if len(accepted_values) >= n_sample:
            break
    return accepted_values, acceptance_rate, i


def metropolis_hastings(likelihood, prior, gen_instr, n_iter=100000, n_sample=100):
    """
    :param likelihood: fct theta -> log( f(x|theta) )
    :param prior: fct theta -> pi(theta)
    :param gen_instr: fct theta -> distr( . |theta)
    """

    thetas = []
    accepted_rate = []

    theta = prior.rvs()
    for i in tqdm(range(n_iter)):
        sampling_instr = gen_instr(theta)  # build the instrumental rho(theta | theta_k)
        candidate = sampling_instr.rvs()  # sample the instr

        candidate_instr = gen_instr(theta)  # build rho(theta_k | candidate)

        log_alpha = likelihood(candidate) + prior.logpdf(candidate) - likelihood(theta) - prior.logpdf(theta) +  sampling_instr.logpdf(candidate) - candidate_instr.logpdf(theta)

        if np.log(np.random.uniform()) < log_alpha:
            theta = candidate
            accepted_rate.append(1)
        else:
            accepted_rate.append(0)
        thetas.append(theta)
        if(sum(accepted_rate)>=n_sample):
            break
    
    return gen_instr(theta), thetas, accepted_rate


def gibbs(params, gen_conditional_posteriors, n_iter=1000):
    """
    :param params: list of initial parameters (theta1_0, ..., thetad_0)
    :param gen_conditional posteriors: list of fct f_i params_k -> distr(theta_i | params_k)  
    """

    samples = [params]
    for i in tqdm(range(n_iter)):
        new_p = []
        for j in range(len(params)):
            thetas = new_p + samples[-1][j:]
            post = gen_conditional_posteriors[j](samples[-1])  # generate the posterior with previous iter
            new_p.append(post.rvs())  # sample the posterior for the next iter
        samples.append(new_p)

    return [gen_conditional_posterior_j(samples[-1]) for gen_conditional_posterior_j in gen_conditional_posteriors], np.array(samples)

def metropolis_hastings_gibbs(likelihood, prior, gen_instr, n_iter=1000):
    """
    :param likelihood: fct theta -> log( f(x|theta) )
    :param prior: fct theta -> pi(theta)
    :param gen_instr: fct theta -> distr( . |theta)
    """

    thetas = list(prior.rvs())
    samples = [thetas]
    for i in range(n_iter):
        new_p = []
        for j in range (len(thetas)):
            thetas = new_p + samples[-1][j:]  # previous theta_ks

            sampling_instr = gen_instr(thetas)  # build the instrumental rho(theta | theta_k)
            candidate = list(sampling_instr.rvs())[j]  # sample the instr

            candidate = new_p + [candidate] + samples[-1][j+1:]  # thetas with the new candidate
            candidate_instr = gen_instr(candidate)  # build rho(theta_k | candidate)

            log_alpha = likelihood(candidate) + prior.logpdf(candidate) - likelihood(thetas) - prior.logpdf(thetas) - sampling_instr.logpdf(candidate) + candidate_instr.logpdf(thetas)

            if np.log(np.random.uniform()) < log_alpha:
                thetas = candidate

            new_p.append(thetas[j])

        samples.append(new_p)
    return np.array(samples)

def def_likelihood(data):
    def f_likelihood(theta):
        return st.weibull_min(c=1)
    return f_likelihood

def def_prior():
    return st.invgamma(a=1)

def def_instr():
    return st.norm()

def def_gen_instr():
    def gen_instr(theta):
        return st.norm(loc=theta, scale=1)
    return gen_instr

if __name__ == '__main__':

    todo = ['ar']
    #todo = ['mh']
    #todo = ['gibbs']

    #  ACCEPTATION REJET
    if 'ar' in todo:
        print('Running Rejection Sampling algorithm...')
        data = st.norm(loc=0).rvs(size=15)
        likelihood = def_likelihood(data)
        prior = def_prior()
        instr = def_instr()

        c = likelihood(data.mean())
        n_sample = 1000
        acc, rate, iter = rejection_sampling(likelihood, prior, instr, c, n_sample=n_sample, n_iter=15000)
        if len(acc)!=n_sample:
            print("Not enough iterations, "+str(len(acc))+" samples generated.")
        print("final iteration = "+str(iter))

        fig, axs = plt.subplots(1,2, figsize=(10,5))
        xx = np.linspace(min(data),max(data),100)
        #axs[0].hist(acc, bins=30, density=True)
        axs[0].set_title('Histogram of sampled posterior distribution')
        sns.distplot(acc,ax = axs[0], label="histogram") #histogramme
        axs[0].plot(xx,instr.pdf(xx),'--',label='instrumental density')
        axs[0].plot(xx,prior.pdf(xx),'--',label='a priori density')
        axs[0].plot(data, np.zeros(len(data)), 'o', label="data")
        axs[0].set_ylim(bottom = -0.2)
        axs[0].legend()

        axs[1].plot(np.cumsum(rate)/range(1, len(rate)+1))
        axs[1].set_title('Acceptance Rate')
        plt.show()
        plt.clf()

    # METROPOLIS HASTINGS
    if 'mh' in todo:
        print('Running Metropolis Hastings algorithm...')
        data = st.norm(loc=0).rvs(size=15)
        likelihood = def_likelihood(data)
        prior = def_prior()
        gen_instr = def_gen_instr()

        posterior, thetas, rate = metropolis_hastings(likelihood, prior, gen_instr)

        fig, axs = plt.subplots(1,3, figsize=(15,5))
        xx = np.linspace(min(data),max(data),100)
        axs[0].set_title('Histogram of sampled posterior distribution')
        sns.distplot(thetas,ax = axs[0], label="histogram") #histogramme
        axs[0].plot(xx,gen_instr(0.5).pdf(xx),'--',label='instrumental density')
        axs[0].plot(xx,prior.pdf(xx),'--',label='a priori density')
        axs[0].plot(data, np.zeros(len(data)), 'o', label="data")
        axs[0].set_ylim(bottom = -0.2)
        axs[0].legend()
        axs[1].plot(range(len(thetas)), thetas)
        axs[1].set_title('Evolution of theta_k')
        axs[2].plot(np.cumsum(rate)/range(1, len(rate)+1))
        axs[2].set_title('Acceptance Rate')
        plt.show()
        plt.clf()

    # GIBBS
    if 'gibbs' in todo:
        print('Running Gibbs algorithm...')
        # exemple du cours, x1.. xn ~ norm(mu, sigma), pi(mu, sigma2) *= 1/sigma2
        
        data = st.norm(loc=5, scale=0.5).rvs(size=1000)
        mean = data.mean()
        n = len(data)

        # on cherhce mu, sigma le calcul donne les posterior conditionnelles
        params = [3, 1]  # mu_0, sigma_0
        posterior_mu = lambda p: st.norm(loc=mean, scale=p[1] / 2)
        posterior_sg = lambda p: st.invgamma(n/2, 0, sum([(xi - p[0])**2 for xi in data])/2)

        posteriors, samples = gibbs(params, [posterior_mu, posterior_sg])
        fig, axs = plt.subplots(1,2, figsize=(10,5))

        axs[0].plot(np.linspace(0, 10, 1000), [posteriors[0].pdf(x) for x in np.linspace(0, 10, 1000)], 'r--')
        axs[0].hist(samples[:, 0], bins=30, density=True)
        axs[0].set_title("Mus")
        axs[1].plot(np.linspace(0, 10, 1000), [posteriors[1].pdf(x) for x in np.linspace(0, 10, 1000)], 'r--')
        axs[1].hist(samples[:, 1], bins=30, density=True)
        axs[1].set_title("Sigmas squared")

    # METROPOLIS HASTINGS WITH GIBBS
    if 'mhwg' in todo:
        print('Running Metropolis Hastings with Gibbs...')
        data = st.norm(loc=5, scale=0.5).rvs(size=10)
        mean = data.mean()
        n = len(data)

        likelihood = lambda p: sum(st.norm(loc=p[0], scale=p[1]).logpdf(data))
        prior = st.multivariate_normal(mean=[8, 5], cov = 1 * np.eye(2))
        instr = lambda p: st.multivariate_normal(mean=p, cov = 1/2 * np.eye(len(p)))

        samples = metropolis_hastings_gibbs(likelihood, prior, instr)

        fig, axs = plt.subplots(1,2, figsize=(10,5))

        #axs[0].plot(np.linspace(0, 10, 1000), [posteriors[0].pdf(x) for x in np.linspace(0, 10, 1000)], 'r--')
        axs[0].hist(samples[:, 0], bins=30, density=True)
        axs[0].set_title("Mus")
        #axs[1].plot(np.linspace(0, 10, 1000), [posteriors[1].pdf(x) for x in np.linspace(0, 10, 1000)], 'r--')
        axs[1].hist(samples[:, 1], bins=30, density=True)
        axs[1].set_title("Sigmas squared")
        plt.show()

