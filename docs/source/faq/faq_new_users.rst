.. _docs_faq_new_users:

.. spelling::

   computecanada
   Niagra
   GMail

New Users
---------

To create a `new CCDB account`_, what should I fill in the form field: 'position'?
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. image:: /img/faq/ccdb-position-formfield.png 

**Response**

Use the following option in the form:

:: 

  external collaborator

For the CCRI field, use the following as input:

:: 

  bws-221-01

What email ID should I use if I am an external collaborator?
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

**Response**

GMail or your work ID should be fine as long as you provide the name of your institution (or college in case of students). 

My account is activated. How do I learn more about Compute Canada Servers and resources available?
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

**Response**

See `Compute Canada Documentation <https://docs.computecanada.ca/wiki/Compute_Canada_Documentation>`_.

I was trying to use `GenPipes deployed in a Docker Container`_. What <tag> value should I use?
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

**Response**

You can use the "latest" tag or one of the tags listed at `GenPipes Docker Hub: <https://hub.docker.com/r/c3genomics/genpipes/tags>`_. If you omit the <tag> Docker will use "latest" by default.

My account is activated but I cannot login into beluga-computecanada or any other node - Cedar, Niagra? What is wrong?
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

**Response**

Check out the current CCDB server status `here <http://status.computecanada.ca/>`_. Many a times, not being able to log in might just be due to system unavailability.

Please note that if you try to log in 3 or more times consecutively with a wrong password, your account gets deactivated and your IP address might get blacklisted. You would need to write to `Compute Canada Support`_ to get that reversed. 

.. _new CCDB account: https://ccdb.computecanada.ca/account_application
.. _GenPipes deployed in a Docker Container: https://genpipes.readthedocs.io/en/latest/deploy/dep_gp_container.html
.. _Compute Canada Support: mailto:support@computecanada.ca