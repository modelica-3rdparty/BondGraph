within ;
package BondGraph

  model Readme
  extends Modelica.Icons.Information;
  annotation (Documentation(info="<html>
  <p>BondGraph is a free Bond Graph library for Modelica that supports 
  graphical representation in Dymola 7.<br><br>
  Copyright (C) 2013 Ilja Alkov, Robin Diekmann
  <br><br>
  BondGraph is free software: you can redistribute it and/or modify it 
  under the terms of the GNU Lesser General Public License as published by the 
  Free Software Foundation, either version 3 of the License, or (at your 
  option) any later version.
  <br><br>
  This program is distributed in the hope that it will be useful, but 
  WITHOUT ANY WARRANTY; without even the implied warranty of 
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser 
  General Public License for more details.
  <br><br>
  You should have received a copy of the GNU Lesser General Public License 
  along with this program. If not, see http://www.gnu.org/licenses.
  <br><br>
  <br><br>
  This software is developed an maintained by<br>
  Ilja Alkov and Robin Diekmann<br>
  mailto: ilja.alkov@fh-bielefeld.de, robin.diekmann@fh-bielefeld.de
  <br><br>
  University of Applied Sciences Bielefeld<br>
  Institute of System Dynamics and Mechatronics<br>
  Bielefeld, Germany<br>
  www.fh-bielefeld.de</p></html>"));
  end Readme;

  package BG_linear

    package Elements

      model C "compliance"

        extends BondGraph.BG_linear.Partial_models.Port_x1;

        parameter Real c = 1 "compliance parameter";
        parameter Real e_0 = 0 "initial effort";

        Real e(start=e_0);
        Real f;

      equation
        port.e = e;
        port.f = f;

        f = der(c*e);

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={125,190,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="C")}));
      end C;

      model Gy "gyrator"

        extends BondGraph.BG_linear.Partial_models.Port_x2;

        parameter Real r = 1 "transformation parameter";

      equation
        port_p.e-r*port_n.f = 0;
        port_n.e-r*port_p.f = 0;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={125,190,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="GY")}));
      end Gy;

      model I "inductance"

        extends BondGraph.BG_linear.Partial_models.Port_x1;

        parameter Real i = 1 "inductance parameter";
        parameter Real f_0 = 0 "initial flow";

        Real e;
        Real f(start=f_0);

      equation
        port.e = e;
        port.f = f;

        e = der(i*f);

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={125,190,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="I")}));
      end I;

      model Jone

        BondGraph.Interfaces.Port_p port_p[p]
          annotation (
            Placement(
              transformation(
                extent={{-100,-10},{-80,10}})));

        BondGraph.Interfaces.Port_n port_n[n]
          annotation (
            Placement(
              transformation(
                extent={{80,-10},{100,10}})));

        parameter Integer p = 0 "number of power inputs" annotation(Dialog(connectorSizing=true));
        parameter Integer n = 0 "number of power outputs" annotation(Dialog(connectorSizing=true));

        Real e[1,p+n];
        Real f[1,p+n];

      equation
        e = cat(2,{port_p.e},-{port_n.e});
        f = cat(2,{port_p.f},{port_n.f});

        if p + n == 0 then
          // do nothing
        else
          zeros(1) = e*ones(p+n);
          zeros(1,p+n-1) = f*cat(1,-ones(1,p+n-1),identity(p+n-1));
        end if;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={125,190,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="1")}));
      end Jone;

      model Jzero

        BondGraph.Interfaces.Port_p port_p[p]
          annotation (
            Placement(
              transformation(
                extent={{-100,-10},{-80,10}})));

        BondGraph.Interfaces.Port_n port_n[n]
          annotation (
            Placement(
              transformation(
                extent={{80,-10},{100,10}})));

        parameter Integer p = 0 "number of power inputs" annotation(Dialog(connectorSizing=true));
        parameter Integer n = 0 "number of power outputs" annotation(Dialog(connectorSizing=true));

        Real e[1,p+n];
        Real f[1,p+n];

      equation
        e = cat(2,{port_p.e},{port_n.e});
        f = cat(2,{port_p.f},-{port_n.f});

        if p + n == 0 then
          // do nothing
        else
          zeros(1) = f*ones(p+n);
          zeros(1,p+n-1) = e*cat(1,-ones(1,p+n-1),identity(p+n-1));
        end if;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={125,190,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="0")}));
      end Jzero;

      model Mc "modulated compliance"

        extends BondGraph.BG_linear.Partial_models.Port_x1;
        extends BondGraph.BG_linear.Partial_models.Port_in_x1;

        Real c "compliance parameter";
        parameter Real e_0 = 0 "initial effort";

        Real e(start=e_0);
        Real f;

      equation
        port.e = e;
        port.f = f;
        port_in = c;

        f = der(c*e);

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={125,190,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="MC")}));
      end Mc;

      model Mgy "modulated gyrator"

        extends BondGraph.BG_linear.Partial_models.Port_x2;
        extends BondGraph.BG_linear.Partial_models.Port_in_x1;

        Real r "transformation parameter";

      equation
        port_p.e-r*port_n.f = 0;
        port_n.e-r*port_p.f = 0;
        port_in = r;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={125,190,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="MGY")}));
      end Mgy;

      model Mi "modulated inductance"

        extends BondGraph.BG_linear.Partial_models.Port_x1;
        extends BondGraph.BG_linear.Partial_models.Port_in_x1;

        Real i "inductance parameter";
        parameter Real f_0 = 0 "initial flow";

        Real e;
        Real f(start=f_0);

      equation
        port.e = e;
        port.f = f;
        port_in = i;

        e = der(i*f);

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={125,190,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="MI")}));
      end Mi;

      model Mr "modulated resistance"

        extends BondGraph.BG_linear.Partial_models.Port_x1;
        extends BondGraph.BG_linear.Partial_models.Port_in_x1;

        Real r "resistance parameter";

      equation
        port.e = r*port.f;
        port_in = r;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={125,190,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="MR")}));
      end Mr;

      model Mse "modulated effort source"

        extends BondGraph.BG_linear.Partial_models.Port_x1;
        extends BondGraph.BG_linear.Partial_models.Port_in_x1;

        Real e_s "source effort";

      equation
        port.e = e_s;
        port_in = e_s;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={125,190,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="MSe")}));
      end Mse;

      model Msf "modulated flow source"

        extends BondGraph.BG_linear.Partial_models.Port_x1;
        extends BondGraph.BG_linear.Partial_models.Port_in_x1;

        Real f_s "source flow";

      equation
        port.f = f_s;
        port_in = f_s;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={125,190,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="MSf")}));
      end Msf;

      model Mtf "modulated transformer"

        extends BondGraph.BG_linear.Partial_models.Port_x2;
        extends BondGraph.BG_linear.Partial_models.Port_in_x1;

        Real r "transformation parameter";

      equation
        port_p.e-r*port_n.e = 0;
        port_n.f-r*port_p.f = 0;
        port_in = r;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={125,190,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="MTF")}));
      end Mtf;

      model R "resistance"

        extends BondGraph.BG_linear.Partial_models.Port_x1;

        parameter Real r = 1 "resistance parameter";

      equation
        port.e = r*port.f;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={125,190,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="R")}));
      end R;

      model Se "effort source"

        extends BondGraph.BG_linear.Partial_models.Port_x1;

        parameter Real e_s = 0 "source effort";

      equation
        port.e = e_s;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={125,190,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="Se")}));
      end Se;

      model Sf "flow source"

        extends BondGraph.BG_linear.Partial_models.Port_x1;

        parameter Real f_s = 0 "source flow";

      equation
        port.f = f_s;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={125,190,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="Sf")}));
      end Sf;

      model Tf "transformer"

        extends BondGraph.BG_linear.Partial_models.Port_x2;

        parameter Real r = 1 "transformation parameter";

      equation
        port_p.e-r*port_n.e = 0;
        port_n.f-r*port_p.f = 0;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={125,190,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="TF")}));
      end Tf;
    end Elements;

    package Sensors

      model Sensor_d_e

        extends BondGraph.BG_linear.Partial_models.Sensor_interface;

        Real y(stateSelect=StateSelect.avoid);

      equation
        port_out = gain*y;

        y = der(e);

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={255,255,0},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="der(e)")}));
      end Sensor_d_e;

      model Sensor_d_f

        extends BondGraph.BG_linear.Partial_models.Sensor_interface;

        Real y(stateSelect=StateSelect.avoid);

      equation
        port_out = gain*y;

        y = der(f);

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={255,255,0},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="der(f)")}));
      end Sensor_d_f;

      model Sensor_d_power

        extends BondGraph.BG_linear.Partial_models.Sensor_interface;

        Real y(stateSelect=StateSelect.avoid);

      equation
        port_out = gain*y;

        y = der(e*f);

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={255,255,0},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="der(P)")}));
      end Sensor_d_power;

      model Sensor_i_e

        extends BondGraph.BG_linear.Partial_models.Sensor_interface;

        parameter Real y_0 = 0 "initial integral value";

        Real y(
          start=y_0,
          stateSelect=StateSelect.prefer);

      equation
        port_out = gain*y;

        der(y) = e;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={255,255,0},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="int(e)")}));
      end Sensor_i_e;

      model Sensor_i_f

        extends BondGraph.BG_linear.Partial_models.Sensor_interface;

        parameter Real y_0 = 0 "initial integral value";

        Real y(
          start=y_0,
          stateSelect=StateSelect.prefer);

      equation
        port_out = gain*y;

        der(y) = f;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={255,255,0},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="int(f)")}));
      end Sensor_i_f;

      model Sensor_i_power

        extends BondGraph.BG_linear.Partial_models.Sensor_interface;

        parameter Real y_0 = 0 "initial integral value";

        Real y(
          start=y_0,
          stateSelect=StateSelect.prefer);

      equation
        port_out = gain*y;

        der(y) = e*f;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={255,255,0},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="int(P)")}));
      end Sensor_i_power;

      model Sensor_p_e

        extends BondGraph.BG_linear.Partial_models.Sensor_interface;

        Real y;

      equation
        port_out = gain*y;

        y = e;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={255,255,0},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="e")}));
      end Sensor_p_e;

      model Sensor_p_f

        extends BondGraph.BG_linear.Partial_models.Sensor_interface;

        Real y;

      equation
        port_out = gain*y;

        y = f;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={255,255,0},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="f")}));
      end Sensor_p_f;

      model Sensor_p_power

        extends BondGraph.BG_linear.Partial_models.Sensor_interface;

        Real y;

      equation
        port_out = gain*y;

        y = e*f;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={255,255,0},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="P")}));
      end Sensor_p_power;
    end Sensors;

    package Partial_models

      partial model Port_in_x1

        BondGraph.Interfaces.Port_in port_in
          annotation (
            Placement(
              transformation(
                extent={{-10,80},{10,100}})));

      end Port_in_x1;

      partial model Port_out_x1

        BondGraph.Interfaces.Port_out port_out
          annotation (
            Placement(
              transformation(
                extent={{-10,-100},{10,-80}})));

      end Port_out_x1;

      partial model Port_x1

        BondGraph.Interfaces.Port port
          annotation (
            Placement(
              transformation(
                extent={{80,-10},{100,10}})));

      end Port_x1;

      partial model Port_x2

        BondGraph.Interfaces.Port_p port_p
          annotation (
            Placement(
              transformation(
                extent={{-100,-10},{-80,10}})));

        BondGraph.Interfaces.Port_n port_n
          annotation (
            Placement(
              transformation(
                extent={{80,-10},{100,10}})));

      end Port_x2;

      partial model Sensor_interface

        extends BondGraph.BG_linear.Partial_models.Port_x2;
        extends BondGraph.BG_linear.Partial_models.Port_out_x1;

        parameter Real gain = 1 "output gain";

        Real e(stateSelect=StateSelect.avoid);
        Real f(stateSelect=StateSelect.avoid);

      equation
        port_p.e-port_n.e = 0;
        port_p.f-port_n.f = 0;
        port_p.e = e;
        port_p.f = f;

      end Sensor_interface;
    end Partial_models;
  end BG_linear;

  package BG_nonlinear
    package Hydraulics

      model Hc "hydraulic capacitance"

        extends BondGraph.BG_linear.Partial_models.Port_x2;

        replaceable BondGraph.BG_nonlinear.Media.Hlp_iso_vg_32.Hprop hprop
          "properties of fluid, model of viscosity and density";
        parameter Real temp = 313.15 "temperature";
        parameter Real part_mass_air = 1.393e-6
          "mass proportion of undissolved air to total mass of mixture";

        parameter Real e_nom = 1e6 "nominal effort";
        parameter Real f_nom = 1e-3 "nominal flow";
        parameter Real c = 1e-8 "compressibility (dv/dp)/v=c=const";
        parameter Real v_0 = 1e-6 "volume at initial pressure";
        parameter Real e_0 = 1e5 "initial pressure";

        Real e(
          nominal=e_nom,
          start=e_0) "effort";
        Real f(nominal=f_nom) "flow";
        Real eta "viscosity";
        Real rho "density";
        Real v(nominal=v_0) "volume of hydraulic node";

      equation
        port_p.e-port_n.e = 0;
        port_p.f-port_n.f = f;
        port_p.e = e;
        hprop.port_in[1] = 0.5*(port_p.e+port_n.e);
        hprop.port_in[2] = temp;
        hprop.port_in[3] = part_mass_air;
        hprop.port_out[1] = eta;
        hprop.port_out[2] = rho;

        v = v_0*exp(c*(e-e_0));
        f = (1/rho)*der(rho*v);

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={181,117,52},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="HC")}));
      end Hc;

      model Hi "hydraulic inductance"

        extends BondGraph.BG_linear.Partial_models.Port_x2;

        replaceable BondGraph.BG_nonlinear.Media.Hlp_iso_vg_32.Hprop hprop
          "properties of fluid, model of viscosity and density";
        parameter Real temp = 313.15 "temperature";
        parameter Real part_mass_air = 1.393e-6
          "mass proportion of undissolved air to total mass of mixture";

        parameter Real e_nom = 1e6 "nominal effort";
        parameter Real f_nom = 1e-3 "nominal flow";
        parameter Real a = 3.14e-4 "cross sectional area of pipe";
        parameter Real l = 0.1 "length of pipe";
        parameter Real f_0 = 0 "initial flow";

        Real e(nominal=e_nom) "effort";
        Real f(
          nominal=f_nom,
          start=f_0) "flow";
        Real eta "viscosity";
        Real rho "density";

      equation
        port_p.e-port_n.e = e;
        port_p.f-port_n.f = 0;
        port_p.f = f;
        hprop.port_in[1] = 0.5*(port_p.e+port_n.e);
        hprop.port_in[2] = temp;
        hprop.port_in[3] = part_mass_air;
        hprop.port_out[1] = eta;
        hprop.port_out[2] = rho;

        e = der((l/a)*rho*f);

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={181,117,52},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="HI")}));
      end Hi;

      model Hr "hydraulic resistance, pipe friction"

        extends BondGraph.BG_linear.Partial_models.Port_x2;

        replaceable BondGraph.BG_nonlinear.Media.Hlp_iso_vg_32.Hprop hprop
          "properties of fluid, model of viscosity and density";
        parameter Real temp = 313.15 "temperature";
        parameter Real part_mass_air = 1.393e-6
          "mass proportion of undissolved air to total mass of mixture";

        parameter Real e_nom = 1e6 "nominal effort";
        parameter Real f_nom = 1e-3 "nominal flow";
        parameter Real a = 3.14e-4 "cross sectional area of pipe";
        parameter Real l = 0.1 "length of pipe";
        parameter Real d_h = 0.01 "hydraulic diameter of pipe d_h";
        parameter Real r_h = 1e-6 "relative hydraulic roughness of pipe k/d_h";
        parameter Real re_min = 1e-3 "minimal reynolds number";
        parameter Real re_crit = 2320 "critical reynolds number";
        parameter Real re_range = 1e1
          "reynolds number range for laminar-turulent transition";
        parameter Real r_min = 1e-2*e_nom/f_nom "minimal resistance";
        parameter Real r_max = 1e6*e_nom/f_nom "maximal resistance";
        parameter Integer par_caus = 2
          "constitutive equation causality, par_caus=1<=>(effort out), par_caus=2<=>(flow out)";

        Real e(nominal=e_nom) "effort";
        Real f(nominal=f_nom) "flow";
        Real r(nominal=e_nom/f_nom) "resistance";
        Real eta "viscosity";
        Real rho "density";
        Real re "reynolds number";
        Real lambda "hydraulic pipe friction factor";

      protected
        parameter Real e_zero = 1e-38*e_nom;
        parameter Real f_zero = 1e-38*f_nom;

      equation
        port_p.e-port_n.e = e;
        port_p.f-port_n.f = 0;
        port_p.f = f;
        hprop.port_in[1] = 0.5*(port_p.e+port_n.e);
        hprop.port_in[2] = temp;
        hprop.port_in[3] = part_mass_air;
        hprop.port_out[1] = eta;
        hprop.port_out[2] = rho;

        e = r*f;

        /* darcy weisbach equation */
        if par_caus == 1 then
          if ((rho*l*lambda)/(2*d_h*a^2))*(abs(f)+f_zero) < r_min then
            r = r_min;
          elseif ((rho*l*lambda)/(2*d_h*a^2))*(abs(f)+f_zero) > r_max then
            r = r_max;
          else
            r = ((rho*l*lambda)/(2*d_h*a^2))*(abs(f)+f_zero);
          end if;
        else
          if (((rho*l*lambda)/(2*d_h*a^2))^(1/2))/((abs(e)+e_zero)^(-1/2)) < r_min then
            r = r_min;
          elseif (((rho*l*lambda)/(2*d_h*a^2))^(1/2))/((abs(e)+e_zero)^(-1/2)) > r_max then
            r = r_max;
          else
            r = (((rho*l*lambda)/(2*d_h*a^2))^(1/2))/((abs(e)+e_zero)^(-1/2));
          end if;
        end if;

        re = (rho*abs(f)*d_h)/(a*eta);

        /* colebrook white equation, alkov diekmann haaland approximation (2012) */
        lambda =
          (64/(re+re_min-re_min*tanh(re/re_min)))*0.5*(-tanh(2.5*(re-re_crit)/re_range)+1)
          +
          ((abs((-1.8/log(10))*log((r_h/3.7)^1.11+6.9/(re+1e-38*re_min)))+1e-38)^(-2))*0.5*(tanh(2.5*(re-re_crit)/re_range)+1);

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={181,117,52},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="HR")}));
      end Hr;

      model Hrl "hydraulic resistance, laminar resistance"

        extends BondGraph.BG_linear.Partial_models.Port_x2;

        replaceable BondGraph.BG_nonlinear.Media.Hlp_iso_vg_32.Hprop hprop
          "properties of fluid, model of viscosity and density";
        parameter Real temp = 313.15 "temperature";
        parameter Real part_mass_air = 1.393e-6
          "mass proportion of undissolved air to total mass of mixture";

        parameter Real e_nom = 1e6 "nominal effort";
        parameter Real f_nom = 1e-3 "nominal flow";
        parameter Real e_ref = 0.5e5 "reference pressure difference";
        parameter Real f_ref = 0.5e-3 "reference volume flow";
        parameter Real eta_ref = 0.0271 "reference viscosity";
        parameter Real p = 1 "volume flow exponent for resistance calculation";
        parameter Real r_min = 1e-2*e_nom/f_nom "minimal resistance";
        parameter Real r_max = 1e6*e_nom/f_nom "maximal resistance";
        parameter Integer par_caus = 2
          "constitutive equation causality, par_caus=1<=>(effort out), par_caus=2<=>(flow out)";

        Real e(nominal=e_nom) "effort";
        Real f(nominal=f_nom) "flow";
        Real r(nominal=e_nom/f_nom) "resistance";
        Real eta "viscosity";
        Real rho "density";

      protected
        parameter Real e_zero = 1e-38*e_nom;
        parameter Real f_zero = 1e-38*f_nom;

      equation
        port_p.e-port_n.e = e;
        port_p.f-port_n.f = 0;
        port_p.f = f;
        hprop.port_in[1] = 0.5*(port_p.e+port_n.e);
        hprop.port_in[2] = temp;
        hprop.port_in[3] = part_mass_air;
        hprop.port_out[1] = eta;
        hprop.port_out[2] = rho;

        e = r*f;

        if par_caus == 1 then
          if ((e_ref*eta)/(eta_ref*abs(f_ref)^p))*((abs(f)+f_zero)^(p-1)) < r_min then
            r = r_min;
          elseif ((e_ref*eta)/(eta_ref*abs(f_ref)^p))*((abs(f)+f_zero)^(p-1)) > r_max then
            r = r_max;
          else
            r = ((e_ref*eta)/(eta_ref*abs(f_ref)^p))*((abs(f)+f_zero)^(p-1));
          end if;
        else
          if (((e_ref*eta)/(eta_ref*abs(f_ref)^p))^(1/p))/((abs(e)+e_zero)^((1/p)-1)) < r_min then
            r = r_min;
          elseif (((e_ref*eta)/(eta_ref*abs(f_ref)^p))^(1/p))/((abs(e)+e_zero)^((1/p)-1)) > r_max then
            r = r_max;
          else
            r = (((e_ref*eta)/(eta_ref*abs(f_ref)^p))^(1/p))/((abs(e)+e_zero)^((1/p)-1));
          end if;
        end if;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={181,117,52},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="HRL")}));
      end Hrl;

      model Hrt "hydraulic resistance, turbulent resistance"

        extends BondGraph.BG_linear.Partial_models.Port_x2;

        replaceable BondGraph.BG_nonlinear.Media.Hlp_iso_vg_32.Hprop hprop
          "properties of fluid, model of viscosity and density";
        parameter Real temp = 313.15 "temperature";
        parameter Real part_mass_air = 1.393e-6
          "mass proportion of undissolved air to total mass of mixture";

        parameter Real e_nom = 1e6 "nominal effort";
        parameter Real f_nom = 1e-3 "nominal flow";
        parameter Real e_ref = 0.5e5 "reference pressure difference";
        parameter Real f_ref = 0.5e-3 "reference volume flow";
        parameter Real rho_ref = 854.307 "reference density";
        parameter Real p = 2 "volume flow exponent for resistance calculation";
        parameter Real r_min = 1e-2*e_nom/f_nom "minimal resistance";
        parameter Real r_max = 1e6*e_nom/f_nom "maximal resistance";
        parameter Integer par_caus = 2
          "constitutive equation causality, par_caus=1<=>(effort out), par_caus=2<=>(flow out)";

        Real e(nominal=e_nom) "effort";
        Real f(nominal=f_nom) "flow";
        Real r(nominal=e_nom/f_nom) "resistance";
        Real eta "viscosity";
        Real rho "density";

      protected
        parameter Real e_zero = 1e-38*e_nom;
        parameter Real f_zero = 1e-38*f_nom;

      equation
        port_p.e-port_n.e = e;
        port_p.f-port_n.f = 0;
        port_p.f = f;
        hprop.port_in[1] = 0.5*(port_p.e+port_n.e);
        hprop.port_in[2] = temp;
        hprop.port_in[3] = part_mass_air;
        hprop.port_out[1] = eta;
        hprop.port_out[2] = rho;

        e = r*f;

        if par_caus == 1 then
          if ((e_ref*rho)/(rho_ref*abs(f_ref)^p))*((abs(f)+f_zero)^(p-1)) < r_min then
            r = r_min;
          elseif ((e_ref*rho)/(rho_ref*abs(f_ref)^p))*((abs(f)+f_zero)^(p-1)) > r_max then
            r = r_max;
          else
            r = ((e_ref*rho)/(rho_ref*abs(f_ref)^p))*((abs(f)+f_zero)^(p-1));
          end if;
        else
          if (((e_ref*rho)/(rho_ref*abs(f_ref)^p))^(1/p))/((abs(e)+e_zero)^((1/p)-1)) < r_min then
            r = r_min;
          elseif (((e_ref*rho)/(rho_ref*abs(f_ref)^p))^(1/p))/((abs(e)+e_zero)^((1/p)-1)) > r_max then
            r = r_max;
          else
            r = (((e_ref*rho)/(rho_ref*abs(f_ref)^p))^(1/p))/((abs(e)+e_zero)^((1/p)-1));
          end if;
        end if;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={181,117,52},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="HRT")}));
      end Hrt;

      model Hse_acc "hydraulic effort source, acceleration pressure"

        extends BondGraph.BG_linear.Partial_models.Port_x2;

        replaceable BondGraph.BG_nonlinear.Media.Hlp_iso_vg_32.Hprop hprop
          "properties of fluid, model of viscosity and density";
        parameter Real temp = 313.15 "temperature";
        parameter Real part_mass_air = 1.393e-6
          "mass proportion of undissolved air to total mass of mixture";

        parameter Real e_nom = 1e6 "nominal effort";
        parameter Real f_nom = 1e-3 "nominal flow";
        parameter Real a = 9.81 "acceleration";
        parameter Real d = 1 "distance in direction of acceleration";

        Real e(nominal=e_nom) "effort";
        Real f(nominal=f_nom) "flow";
        Real eta "viscosity";
        Real rho "density";

      equation
        port_p.e-port_n.e = e;
        port_p.f-port_n.f = 0;
        port_p.f = f;
        hprop.port_in[1] = 0.5*(port_p.e+port_n.e);
        hprop.port_in[2] = temp;
        hprop.port_in[3] = part_mass_air;
        hprop.port_out[1] = eta;
        hprop.port_out[2] = rho;

        e = -rho*a*d; // acceleration pressure

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={181,117,52},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="HSe
acc")}));
      end Hse_acc;

      model Hse_ind
        "hydraulic effort source, inductance pressure, cross section change"

        extends BondGraph.BG_linear.Partial_models.Port_x2;

        replaceable BondGraph.BG_nonlinear.Media.Hlp_iso_vg_32.Hprop hprop
          "properties of fluid, model of viscosity and density";
        parameter Real temp = 313.15 "temperature";
        parameter Real part_mass_air = 1.393e-6
          "mass proportion of undissolved air to total mass of mixture";

        parameter Real e_nom = 1e6 "nominal effort";
        parameter Real f_nom = 1e-3 "nominal flow";
        parameter Real a_1 = 1 "inlet cross section";
        parameter Real a_2 = 1 "outlet cross section";

        Real e(nominal=e_nom) "effort";
        Real f(nominal=f_nom) "flow";
        Real eta "viscosity";
        Real rho "density";

      equation
        port_p.e-port_n.e = e;
        port_p.f-port_n.f = 0;
        port_p.f = f;
        hprop.port_in[1] = 0.5*(port_p.e+port_n.e);
        hprop.port_in[2] = temp;
        hprop.port_in[3] = part_mass_air;
        hprop.port_out[1] = eta;
        hprop.port_out[2] = rho;

        e = -0.5*rho*abs(f)^2*(a_1^(-2)-a_2^(-2)); // bernoulli relation

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={181,117,52},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="HSe
ind")}));
      end Hse_ind;

      model Mhc "modulated hydraulic capacitance"

        extends BondGraph.BG_linear.Partial_models.Port_x2;

        BondGraph.Interfaces.Port_in port_in
          annotation (
            Placement(
              transformation(
                extent={{-10,80},{10,100}})));

        replaceable BondGraph.BG_nonlinear.Media.Hlp_iso_vg_32.Hprop hprop
          "properties of fluid, model of viscosity and density";
        parameter Real temp = 313.15 "temperature";
        parameter Real part_mass_air = 1.393e-6
          "mass proportion of undissolved air to total mass of mixture";

        parameter Real e_nom = 1e6 "nominal effort";
        parameter Real f_nom = 1e-3 "nominal flow";
        parameter Real c = 1e-8 "compressibility (dv/dp)/v";
        Real v_0 "volume at initial pressure";
        parameter Real e_0 = 1e5 "initial pressure";

        Real e(
          nominal=e_nom,
          start=e_0) "effort";
        Real f(nominal=f_nom) "flow";
        Real eta "viscosity";
        Real rho "density";
        Real v "volume of hydraulic node";

      equation
        port_p.e-port_n.e = 0;
        port_p.f-port_n.f = f;
        port_p.e = e;
        hprop.port_in[1] = 0.5*(port_p.e+port_n.e);
        hprop.port_in[2] = temp;
        hprop.port_in[3] = part_mass_air;
        hprop.port_out[1] = eta;
        hprop.port_out[2] = rho;
        port_in = v_0;

        v = v_0*exp(c*(e-e_0));
        f = (1/rho)*der(rho*v);

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={181,117,52},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="MHC")}));
      end Mhc;

      model Mhrt "signal controlled hydraulic resistance, turbulent resistance"

        extends BondGraph.BG_linear.Partial_models.Port_x2;
        extends BondGraph.BG_linear.Partial_models.Port_in_x1;

        replaceable BondGraph.BG_nonlinear.Media.Hlp_iso_vg_32.Hprop hprop
          "properties of fluid, model of viscosity and density";
        parameter Real temp = 313.15 "temperature";
        parameter Real part_mass_air = 1.393e-6
          "mass proportion of undissolved air to total mass of mixture";

        parameter Real e_nom = 1e6 "nominal effort";
        parameter Real f_nom = 1e-3 "nominal flow";
        parameter Real e_ref = 0.5e5 "reference pressure difference";
        parameter Real f_ref = 0.5e-3 "reference volume flow";
        parameter Real rho_ref = 854.307 "reference density";
        parameter Real p = 2 "volume flow exponent for resistance calculation";
        parameter Real r_min = 1e-2*e_nom/f_nom "minimal resistance";
        parameter Real r_max = 1e6*e_nom/f_nom "maximal resistance";
        parameter Integer par_caus = 2
          "constitutive equation causality, par_caus=1<=>(effort out), par_caus=2<=>(flow out)";

        Real e(nominal=e_nom) "effort";
        Real f(nominal=f_nom) "flow";
        Real r(nominal=e_nom/f_nom) "resistance";
        Real eta "viscosity";
        Real rho "density";
        Real s_c(nominal=1) "control signal";

      protected
        parameter Real e_zero = 1e-38*e_nom;
        parameter Real f_zero = 1e-38*f_nom;

      equation
        port_p.e-port_n.e = e;
        port_p.f-port_n.f = 0;
        port_p.f = f;
        hprop.port_in[1] = 0.5*(port_p.e+port_n.e);
        hprop.port_in[2] = temp;
        hprop.port_in[3] = part_mass_air;
        hprop.port_out[1] = eta;
        hprop.port_out[2] = rho;
        port_in = s_c;

        e = (1/s_c)*r*f;

        if par_caus == 1 then
          if ((e_ref*rho)/(rho_ref*abs(f_ref)^p))*((abs(f)+f_zero)^(p-1)) < r_min then
            r = r_min;
          elseif ((e_ref*rho)/(rho_ref*abs(f_ref)^p))*((abs(f)+f_zero)^(p-1)) > r_max then
            r = r_max;
          else
            r = ((e_ref*rho)/(rho_ref*abs(f_ref)^p))*((abs(f)+f_zero)^(p-1));
          end if;
        else
          if (((e_ref*rho)/(rho_ref*abs(f_ref)^p))^(1/p))/((abs(e)+e_zero)^((1/p)-1)) < r_min then
            r = r_min;
          elseif (((e_ref*rho)/(rho_ref*abs(f_ref)^p))^(1/p))/((abs(e)+e_zero)^((1/p)-1)) > r_max then
            r = r_max;
          else
            r = (((e_ref*rho)/(rho_ref*abs(f_ref)^p))^(1/p))/((abs(e)+e_zero)^((1/p)-1));
          end if;
        end if;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={181,117,52},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="MHRT")}));
      end Mhrt;
    end Hydraulics;

    package Mechanics

      model Mse_masy "mechanical effort source, asynchronous electric motor"

        BondGraph.Interfaces.Port port
          annotation (
            Placement(
              transformation(
                extent={{80,-10},{100,10}})));

        BondGraph.Interfaces.Port_in port_in
          annotation (
            Placement(
              transformation(
                extent={{-10,80},{10,100}})));

      /* static torque values */
        parameter Real par_torq_1 = 120.78 "static torque, start value";
        parameter Real par_torq_2 = 115 "static torque, saddle value";
        parameter Real par_torq_3 = 120.78 "static torque, breakdown value";
        parameter Real par_torq_4 = 36.6 "static torque, nominal value";

      /* rotational frequency values */
        parameter Real par_freq_1 = 0 "rotational frequency, start value";
        parameter Real par_freq_2 = 5 "rotational frequency, saddle value";
        parameter Real par_freq_3 = 15 "rotational frequency, breakdown value";
        parameter Real par_freq_4 = 1435/60
          "rotational frequency, nominal value";

      /* torque build-up */
        parameter Real par_dyn_1_1 = 0
          "initial torque, relative value, par_dyn_1_1 in R and par_dyn_1_1 in [0,1]";
        parameter Real par_dyn_1_2 = 50e-3 "torque bild-up, time";
        parameter Real par_dyn_1_3 = 0.99
          "torque bild-up, relative value, par_dyn_1_3 in R and par_dyn_1_3 in (0,1)";

      /* torque starting oscillation */
        parameter Real par_dyn_2_1 = par_torq_4
          "starting oscillation, torque amplitude";
        parameter Real par_dyn_2_2 = 0.02
          "starting oscillation, damping parameter, par_dyn_2_2 in R and par_dyn_2_2 in (0,1)";
        parameter Real par_dyn_2_3 = 50 "starting oscillation, frequency";

        Real e "effort";
        Real f "flow";
        Real t_stat "static torque";
        Real f_rot "rotational frequency";
        Real u "input variable, u in R and u in [0,+inf)";
        Real x_dyn[3](start={par_dyn_1_1,0,par_dyn_1_1*par_dyn_2_1*(1-par_dyn_2_2^2)/(4*asin(1)*par_dyn_2_3)})
          "dynamic state";

      equation
        port.e = e;
        port.f = f;
        port_in = u;

        {{-par_dyn_1_2/log(1-par_dyn_1_3), 0, 0},
         {0, 0, 1},
         {0, (1-par_dyn_2_2^2)/(4*asin(1)*par_dyn_2_3)^2, 0}}
        *
        der(x_dyn)
        +
        {{1, 0, 0},
         {0, -1, 0},
         {0, 2*par_dyn_2_2*sqrt(1-par_dyn_2_2^2)/(4*asin(1)*par_dyn_2_3), 1}}
        *
        x_dyn
        =
        {1, 0, par_dyn_2_1*(1-par_dyn_2_2^2)/(4*asin(1)*par_dyn_2_3)}
        *
        u;

        e = t_stat*x_dyn[1] + x_dyn[2];

        f = 4*asin(1)*f_rot;

        if f_rot < par_freq_2 then
          t_stat =
          ((par_torq_1-par_torq_2)/(par_freq_1-par_freq_2)^2)
          *
          ((f_rot - par_freq_2)^2)
          +
          par_torq_2;
        elseif f_rot >= par_freq_2 and f_rot < (par_freq_2+par_freq_3)/2 then
          t_stat =
          (((par_torq_2+par_torq_3)/2-par_torq_2)/((par_freq_2+par_freq_3)/2-par_freq_2)^2)
          *
          ((f_rot-par_freq_2)^2)
          +
          par_torq_2;
        elseif f_rot >= (par_freq_2+par_freq_3)/2 and f_rot < par_freq_3 then
          t_stat =
          (((par_torq_2+par_torq_3)/2-par_torq_3)/((par_freq_2+par_freq_3)/2-par_freq_3)^2)
          *
          ((f_rot-par_freq_3)^2)
          +
          par_torq_3;
        else
          t_stat =
          ((par_torq_4-par_torq_3)/(par_freq_4-par_freq_3)^2)
          *
          ((f_rot-par_freq_3)^2)
          +
          par_torq_3;
        end if;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={181,117,52},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="MSe
masy")}));
      end Mse_masy;

      model Mse_masy_stat
        "mechanical effort source, asynchronous electric motor, stationary characteristic"

        BondGraph.Interfaces.Port port
          annotation (
            Placement(
              transformation(
                extent={{80,-10},{100,10}})));

        BondGraph.Interfaces.Port_in port_in
          annotation (
            Placement(
              transformation(
                extent={{-10,80},{10,100}})));

      /* static torque values */
        parameter Real par_torq_1 = 120.78 "static torque, start value";
        parameter Real par_torq_2 = 115 "static torque, saddle value";
        parameter Real par_torq_3 = 120.78 "static torque, breakdown value";
        parameter Real par_torq_4 = 36.6 "static torque, nominal value";

      /* rotational frequency values */
        parameter Real par_freq_1 = 0 "rotational frequency, start value";
        parameter Real par_freq_2 = 5 "rotational frequency, saddle value";
        parameter Real par_freq_3 = 15 "rotational frequency, breakdown value";
        parameter Real par_freq_4 = 1435/60
          "rotational frequency, nominal value";

        Real e "effort";
        Real f "flow";
        Real t_stat "static torque";
        Real f_rot "rotational frequency";
        Real u "input variable, u in R and u in [0,+inf)";

      equation
        port.e = e;
        port.f = f;
        port_in = u;

        e = t_stat*u;

        f = 4*asin(1)*f_rot;

        if f_rot < par_freq_2 then
          t_stat =
          ((par_torq_1-par_torq_2)/(par_freq_1-par_freq_2)^2)
          *
          ((f_rot - par_freq_2)^2)
          +
          par_torq_2;
        elseif f_rot >= par_freq_2 and f_rot < (par_freq_2+par_freq_3)/2 then
          t_stat =
          (((par_torq_2+par_torq_3)/2-par_torq_2)/((par_freq_2+par_freq_3)/2-par_freq_2)^2)
          *
          ((f_rot-par_freq_2)^2)
          +
          par_torq_2;
        elseif f_rot >= (par_freq_2+par_freq_3)/2 and f_rot < par_freq_3 then
          t_stat =
          (((par_torq_2+par_torq_3)/2-par_torq_3)/((par_freq_2+par_freq_3)/2-par_freq_3)^2)
          *
          ((f_rot-par_freq_3)^2)
          +
          par_torq_3;
        else
          t_stat =
          ((par_torq_4-par_torq_3)/(par_freq_4-par_freq_3)^2)
          *
          ((f_rot-par_freq_3)^2)
          +
          par_torq_3;
        end if;

        annotation (
          Icon(
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={181,117,52},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="MSe
masys")}));
      end Mse_masy_stat;
    end Mechanics;

    package Media

      package Generic

        function eta_air
          "viscosity of air, approximation function of experimental data"

          input Real press;
          input Real temp;
          output Real eta;

        protected
          parameter Real press_ref = 1.0e5 "reference pressure";
          parameter Real temp_ref = 288.15 "reference temperature";
          parameter Real par_eta_press = 1
            "viscosity index for pressure dependency of air";
          parameter Real par_eta_temp = 0.132
            "viscosity index for temperature dependency of air";
          parameter Real eta_ref = 16.726e-6
            "viscosity at reference conditions of air";

        algorithm
          eta := eta_ref
                 *
                 ( ( par_eta_press*(press^0.01)*(temp^0.75)
                     +
                     par_eta_temp*press*(temp^(-2)))
                   /
                   ( par_eta_press*(press_ref^0.01)*(temp_ref^0.75)
                     +
                     par_eta_temp*press_ref*(temp_ref^(-2))));

        end eta_air;

        function eta_mix "viscosity of mixture"

          input Real parts_vol[:];
          input Real eta_parts[size(parts_vol,1)];
          output Real eta;

        algorithm
          eta := eta_parts
                 *
                 ( (parts_vol .^ (2/3))
                   /
                   (ones(size(parts_vol,1)) * (parts_vol .^ (2/3))));

        end eta_mix;

        function eta_roelands
          "roelands relation, viscosity, pressure, temperature"

          input Real press;
          input Real temp;
          input Real press_ref;
          input Real temp_ref;
          input Real par_eta_press;
          input Real par_eta_temp;
          input Real eta_ref;
          output Real eta;

        algorithm
          eta := eta_ref
                 *
                 exp(
                   log(eta_ref/(6.315e-5))
                   *
                   (-1+((1+(press-press_ref)/(1.96e8))^par_eta_press)
                   *
                   (((temp-138)/(temp_ref-138))^(-par_eta_temp))));

        end eta_roelands;

        function parts_vol "volume parts of mass components"

          input Real parts_mass[:];
          input Real rho_parts[size(parts_mass,1)];
          output Real parts[size(parts_mass,1)];

        algorithm
          parts := (parts_mass .* (1 ./ rho_parts))
                   /
                   (parts_mass * (1 ./ rho_parts));

        end parts_vol;

        function press_cav
          "cavitation pressure of liquid fluid, approximation function of experimental data"

          input Real temp;
          input Real temp_ref;
          input Real par_press_cav_temp;
          input Real press_cav_ref;
          output Real press;

        algorithm
          press := press_cav_ref
                   *
                   (temp/temp_ref)
                   *
                   exp(
                     par_press_cav_temp*(temp-temp_ref));

        end press_cav;

        function rho_air "density of air, ideal gas law"

          input Real press;
          input Real temp;
          output Real rho;

        protected
          parameter Real press_ref = 1.0e5 "reference pressure";
          parameter Real temp_ref = 288.15 "reference temperature";
          parameter Real par_rho_press = 1
            "density index for pressure dependency of air";
          parameter Real par_rho_temp = -1
            "density index for temperature dependency of air";
          parameter Real rho_ref = 1.211
            "density at reference conditions of air";

        algorithm
          rho := rho_ref
                 *
                 ((press/press_ref)^par_rho_press)
                 *
                 ((temp/temp_ref)^par_rho_temp);

        end rho_air;

        function rho_exp
          "exponential relation, density, pressure, temperature, compressibility factor, coefficient of thermal expansion"

          input Real press;
          input Real temp;
          input Real press_ref;
          input Real temp_ref;
          input Real par_rho_press;
          input Real par_rho_temp;
          input Real rho_ref;
          output Real rho;

        algorithm
          rho := rho_ref
                 *
                 exp(
                   par_rho_press*(press-press_ref)
                   -
                   par_rho_temp*(temp-temp_ref));

        end rho_exp;

        function rho_mix "density of mixture"

          input Real parts_mass[:];
          input Real rho_parts[size(parts_mass,1)];
          output Real rho;

        algorithm
          rho := 1
                 /
                 (parts_mass * (1 ./ rho_parts));

        end rho_mix;
      end Generic;

      package Hlp_iso_vg_32

        package Hlp_iso_vg_32_par

          constant Real press_ref = 1.0e5 "reference pressure";
          constant Real temp_ref = 288.15 "reference temperature";

          constant Real par_eta_press_oil = 0.6
            "viscosity index for pressure dependency (roelands relation) of oil";
          constant Real par_eta_temp_oil = 1.196
            "viscosity index for temperature dependency (roelands relation) of oil";
          constant Real eta_ref_oil = 93.413e-3
            "viscosity at reference conditions of oil";

          constant Real par_rho_press_oil = 7.142e-10
            "density index for pressure dependency (compressibility factor) of oil";
          constant Real par_rho_temp_oil = 7.0e-4
            "density index for temperature dependency (coefficient of thermal expansion) of oil";
          constant Real rho_ref_oil = 870.0
            "density at reference conditions of oil";

          constant Real par_press_cav_temp_oil = 20.68e-3
            "cavitation pressure index for temperature dependency of oil";
          constant Real press_cav_ref_oil = 88.64
            "cavitation pressure at reference temperature of oil";

          constant Real part_mass_ref_air = 1.393e-6
            "mass proportion of undissolved air to total mass of mixture at reference conditions";

        end Hlp_iso_vg_32_par;

        model Hprop
          "properties of fluid, viscosity and density of oil-air mixture"

          BondGraph.Interfaces.Port_in port_in[3]
            annotation (
              Placement(
                transformation(
                  extent={{-10,80},{10,100}})));

          BondGraph.Interfaces.Port_out port_out[2]
            annotation (
              Placement(
                transformation(
                  extent={{-10,-100},{10,-80}})));

          Real press(start=1e5);
          Real temp;
          Real part_mass_air;
          Real parts_mass[2];
          Real parts_vol[2];
          Real eta_air;
          Real eta_oil;
          Real rho_air;
          Real rho_oil;
          Real eta_mix(nominal=eta_ref_oil);
          Real rho_mix(nominal=rho_ref_oil);

        protected
          constant Real press_ref = 1.0e5 "reference pressure";
          constant Real temp_ref = 288.15 "reference temperature";

          constant Real par_rho_press_air = 1
            "density index for pressure dependency of air";
          constant Real par_rho_temp_air = -1
            "density index for temperature dependency of air";
          constant Real rho_ref_air = 1.211
            "density at reference conditions of air";
          constant Real par_eta_press_air = 1
            "viscosity index for pressure dependency of air";
          constant Real par_eta_temp_air = 0.132
            "viscosity index for temperature dependency of air";
          constant Real eta_ref_air = 16.726e-6
            "viscosity at reference conditions of air";

          constant Real par_eta_press_oil = 0.6
            "viscosity index for pressure dependency (roelands relation) of oil";
          constant Real par_eta_temp_oil = 1.196
            "viscosity index for temperature dependency (roelands relation) of oil";
          constant Real eta_ref_oil = 93.413e-3
            "viscosity at reference conditions of oil";
          constant Real par_rho_press_oil = 7.142e-10
            "density index for pressure dependency (compressibility factor) of oil";
          constant Real par_rho_temp_oil = 7.0e-4
            "density index for temperature dependency (coefficient of thermal expansion) of oil";
          constant Real rho_ref_oil = 870.0
            "density at reference conditions of oil";

          constant Real par_press_cav_temp_oil = 20.68e-3
            "cavitation pressure index for temperature dependency of oil";
          constant Real press_cav_ref_oil = 88.64
            "cavitation pressure at reference temperature of oil";

          constant Real part_mass_ref_air = 1.393e-6
            "mass proportion of undissolved air to total mass of mixture at reference conditions";

        equation
          press = if port_in[1] < 1e0 then 1e0 elseif port_in[1] > 1e8 then 1e8 else port_in[1];
          temp = if port_in[2] < 1e2 then 1e2 elseif port_in[2] > 1e3 then 1e3 else port_in[2];
          part_mass_air = if port_in[3] < 0e0 then 0e0 elseif port_in[3] > 1e0 then 1e0 else port_in[3];
          port_out = {eta_mix,rho_mix};

          parts_mass = {1-part_mass_air,part_mass_air};

          rho_air = rho_ref_air
                    *
                    ((press/press_ref)^par_rho_press_air)
                    *
                    ((temp/temp_ref)^par_rho_temp_air);

          rho_oil = rho_ref_oil
                    *
                    exp(
                      par_rho_press_oil*(press-press_ref)
                      -
                      par_rho_temp_oil*(temp-temp_ref));

          rho_mix = 1/(parts_mass * (1 ./ {rho_oil,rho_air}));

          parts_vol = (parts_mass .* (1 ./ {rho_oil,rho_air}))
                      /
                      (parts_mass * (1 ./ {rho_oil,rho_air}));

          eta_air = eta_ref_air
                    *
                    ((par_eta_press_air*(press^0.01)*(temp^0.75)
                      +
                      par_eta_temp_air*press*(temp^(-2)))
                     /
                     (par_eta_press_air*(press_ref^0.01)*(temp_ref^0.75)
                      +
                      par_eta_temp_air*press_ref*(temp_ref^(-2))));

          eta_oil = eta_ref_oil
                    *
                    exp(
                      log(eta_ref_oil/(6.315e-5))
                      *
                      (-1+((1+(press-press_ref)/(1.96e8))^par_eta_press_oil)
                      *
                      (((temp-138)/(temp_ref-138))^(-par_eta_temp_oil))));

          eta_mix = {eta_oil,eta_air}
                    *
                    ((parts_vol .^ (2/3))
                     /
                     (ones(size(parts_vol,1)) * (parts_vol .^ (2/3))));

          annotation (
            Icon(
              graphics={Ellipse(
                  extent={{-100,100},{100,-100}},
                  fillColor={181,117,52},
                  fillPattern=FillPattern.Solid,
                  lineColor={0,0,0}), Text(
                  extent={{0,30},{0,-30}},
                  pattern=LinePattern.None,
                  fillColor={85,170,255},
                  fillPattern=FillPattern.Solid,
                  lineColor={0,0,0},
                  textString="hprop")}));
        end Hprop;
      end Hlp_iso_vg_32;
    end Media;
  end BG_nonlinear;

  package Examples

    package Valves

      model Switch_2x2

        parameter Real par_c_1[3] = (1/0.1e5)*{-1,1,-1}
          "control gains, {u, p_1, p_2}";
        parameter Real par_c_2 = 0
          "control signal for par_c_1*{u, p_1, p_2} = 0";
        parameter Real par_dyn_1 = 1e6 "control signal gain";
        parameter Real par_dyn_2 = -1e6 "state gain";
        parameter Real par_dyn_3 = 0 "initial state";
        parameter Real par_dyn_4 = -1/10e-3 "minimal der(state) value";
        parameter Real par_dyn_5 = 1/10e-3 "maximal der(state) value";
        parameter Real par_stat_1[2] = {0,1} "positions assignment";
        parameter Real par_stat_2[2] = {1,1}
          "coverage parameter for state < p_x[i]";
        parameter Real par_stat_3[2] = {1,1}
          "coverage parameter for state > p_x[i]";
        parameter Real par_res_1 = 313.15 "temperature";
        parameter Real par_res_2 = 1.393e-6
          "mass proportion of undissolved air to total mass of mixture";
        parameter Real par_res_3 = 1e6 "nominal effort";
        parameter Real par_res_4 = 1e-3 "nominal flow";
        parameter Real par_res_5 = 0.5e5 "reference pressure difference";
        parameter Real par_res_6 = 0.5e-3 "reference volume flow";
        parameter Real par_res_7 = 854.307 "reference density";
        parameter Real par_res_8 = 2
          "volume flow exponent for resistance calculation";
        parameter Real par_res_9 = 1e-2*par_res_3/par_res_4
          "minimal resistance";
        parameter Real par_res_10 = 1e6*par_res_3/par_res_4
          "maximal resistance";
        parameter Integer par_res_11 = 2
          "constitutive equation causality, par_res_5=1(effort out), par_res_5=2(flow out)";

        BondGraph.Interfaces.Port_p port_p
          annotation (
            Placement(
              transformation(
                extent={{-200,-10},{-180,10}}),
              iconTransformation(
                extent={{-100,-10},{-80,10}})));

        BondGraph.Interfaces.Port_n port_n
          annotation (
            Placement(
              transformation(
                extent={{180,-10},{200,10}}),
              iconTransformation(
                extent={{80,-10},{100,10}})));

        BondGraph.Interfaces.Port_in port_in
          annotation (
            Placement(
              transformation(
                extent={{-10,180},{10,200}}),
              iconTransformation(
                extent={{-10,180},{10,200}})));

        BondGraph.BG_linear.Sensors.Sensor_p_e sensor_p_e1
          annotation (
            Placement(
              transformation(
                extent={{-160,-10},{-140,10}})));

        BondGraph.BG_linear.Sensors.Sensor_p_e sensor_p_e2
          annotation (
            Placement(
              transformation(
                extent={{140,-10},{160,10}})));

        BondGraph.BG_linear.Elements.Jzero jzero1(p=1, n=1)
          annotation (
            Placement(
              transformation(
                extent={{-120,-10},{-100,10}})));

        BondGraph.BG_linear.Elements.Jzero jzero2(p=1, n=1)
          annotation (
            Placement(
              transformation(
                extent={{100,-10},{120,10}})));

        BondGraph.Blocks.Const const1(c={{par_c_2}})
          annotation (
            Placement(
              transformation(
                extent={{-50,180},{-30,200}})));

        BondGraph.Blocks.Prop_sat prop_sat1(
          b={cat(   1,
                    {1},
                    par_c_1)},
          x_min={par_stat_1[1]},
          x_max={par_stat_1[size(par_stat_1, 1)]})
          annotation (
            Placement(
              transformation(
                extent={{-10,120},{10,140}})));

        BondGraph.Blocks.Prop_sat prop_sat2(
          b=identity(size(par_stat_1, 1)),
          x_min=1e-9*ones(size(par_stat_1, 1)),
          x_max=ones(size(par_stat_1, 1)))
          annotation (
            Placement(
              transformation(
                extent={{-10,30},{10,50}})));

        BondGraph.Blocks.Lin_sat lin_sat1(
          b={{par_dyn_1}},
          a={{par_dyn_2}},
          x_0={par_dyn_3},
          x_min={par_stat_1[1]},
          x_max={par_stat_1[size(par_stat_1, 1)]},
          x_der_min={par_dyn_4},
          x_der_max={par_dyn_5})
          annotation (
            Placement(
              transformation(
                extent={{-10,90},{10,110}})));

        BondGraph.Blocks.Sig_c sig_c1(
          p_x=par_stat_1,
          p_cn=par_stat_2,
          p_cp=par_stat_3)
          annotation (
            Placement(
              transformation(
                extent={{-10,60},{10,80}})));

        BondGraph.BG_nonlinear.Hydraulics.MHRT mhrt1(
          temp=par_res_1,
          part_mass_air=par_res_2,
          e_nom=par_res_3,
          f_nom=par_res_4,
          e_ref=par_res_5,
          f_ref=par_res_6,
          rho_ref=par_res_7,
          p=par_res_8,
          r_min=par_res_9,
          r_max=par_res_10,
          par_caus=par_res_11)
          annotation (
            Placement(
              transformation(
                extent={{-10,-10},{10,10}})));

      equation
        connect(const1.port_out[1, 1], prop_sat1.port_in[1]) annotation (Line(
            points={{-40,181},{-40,170},{0,170},{0,139}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(port_in, prop_sat1.port_in[2]) annotation (Line(
            points={{0,190},{0,139}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(prop_sat1.port_out, lin_sat1.port_in) annotation (Line(
            points={{0,121},{0,109}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(lin_sat1.port_out[1], sig_c1.port_in) annotation (Line(
            points={{0,91},{0,79}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sig_c1.port_out, prop_sat2.port_in) annotation (Line(
            points={{0,61},{0,49}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(port_p, sensor_p_e1.port_p) annotation (Line(
            points={{-190,0},{-159,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_p_e2.port_n, port_n) annotation (Line(
            points={{159,0},{190,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_p_e1.port_n, jzero1.port_p[1]) annotation (Line(
            points={{-141,0},{-119,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero2.port_n[1], sensor_p_e2.port_p) annotation (Line(
            points={{119,0},{141,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_p_e1.port_out, prop_sat1.port_in[3]) annotation (Line(
            points={{-150,-9},{-150,-20},{-130,-20},{-130,170},{0,170},{0,139}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_p_e2.port_out, prop_sat1.port_in[4]) annotation (Line(
            points={{150,-9},{150,-20},{130,-20},{130,170},{0,170},{0,139}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(prop_sat2.port_out[2], mhrt1.port_in) annotation (Line(
            points={{0,31},{0,9}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero1.port_n[1], mhrt1.port_p) annotation (Line(
            points={{-101,0},{-9,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(mhrt1.port_n, jzero2.port_p[1]) annotation (Line(
            points={{9,0},{101,0}},
            color={0,0,0},
            smooth=Smooth.None));

        annotation (
          Diagram(
            coordinateSystem(
              preserveAspectRatio=true,
              extent={{-200,-200},{200,200}}),
            graphics),
          Icon(
            coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-200},{100,200}}),
            graphics={
                Rectangle(
                extent={{-100,200},{100,-200}},
                pattern=LinePattern.None,
                fillColor={0,0,0},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}),
                Rectangle(
                extent={{-80,180},{80,-180}},
                pattern=LinePattern.None,
                fillColor={181,117,52},
                fillPattern=FillPattern.Solid),
                        Text(
                  extent={{0,30},{0,-30}},
                  lineColor={0,0,0},
                  lineThickness=1,
                  fillColor={51,249,65},
                  fillPattern=FillPattern.Solid,
                textString="switch
2x2
")}));
      end Switch_2x2;

      model Switch_3x2

        parameter Real par_c_1[4] = (1/0.1e5)*{-1,0,0,1}
          "control gains, {u, p_1, p_2, p_3}";
        parameter Real par_c_2 = 0
          "control signal for par_c_1*{u, p_1, p_2, p_3} = 0";
        parameter Real par_dyn_1 = 1e6 "control signal gain";
        parameter Real par_dyn_2 = -1e6 "state gain";
        parameter Real par_dyn_3 = 0 "initial state";
        parameter Real par_dyn_4 = -1/10e-3 "minimal der(state) value";
        parameter Real par_dyn_5 = 1/10e-3 "maximal der(state) value";
        parameter Real par_stat_1[2] = {0,1} "positions assignment";
        parameter Real par_stat_2[2] = {0.5,0.5}
          "coverage parameter for state < p_x[i]";
        parameter Real par_stat_3[2] = {0.5,0.5}
          "coverage parameter for state > p_x[i]";
        parameter Real par_res_1[2] = 313.15*{1,1} "temperature";
        parameter Real par_res_2[2] = 1.393e-6*{1,1}
          "mass proportion of undissolved air to total mass of mixture";
        parameter Real par_res_3[2] = 1e6*{1,1} "nominal effort";
        parameter Real par_res_4[2] = 1e-3*{1,1} "nominal flow";
        parameter Real par_res_5[2] = 0.5e5*{1,1}
          "reference pressure difference";
        parameter Real par_res_6[2] = 0.5e-3*{1,1} "reference volume flow";
        parameter Real par_res_7[2] = 854.307*{1,1} "reference density";
        parameter Real par_res_8[2] = 2*{1,1}
          "volume flow exponent for resistance calculation";
        parameter Real par_res_9[2] = 1e-2*par_res_3 ./ par_res_4
          "minimal resistance";
        parameter Real par_res_10[2] = 1e6*par_res_3 ./ par_res_4
          "maximal resistance";
        parameter Integer par_res_11[2] = 2*{1,1}
          "constitutive equation causality, par_res_5=1(effort out), par_res_5=2(flow out)";

        BondGraph.Interfaces.Port_p port_p1
          annotation (
            Placement(
              transformation(
                extent={{-200,-10},{-180,10}}),
              iconTransformation(
                extent={{-100,90},{-80,110}})));

        BondGraph.Interfaces.Port_p port_p2
          annotation (
            Placement(
              transformation(
                extent={{-200,-50},{-180,-30}}),
              iconTransformation(
                extent={{-100,-110},{-80,-90}})));

        BondGraph.Interfaces.Port_n port_n
          annotation (
            Placement(
              transformation(
                extent={{180,-10},{200,10}}),
              iconTransformation(
                extent={{80,-10},{100,10}})));

        BondGraph.Interfaces.Port_in port_in
          annotation (
            Placement(
              transformation(
                extent={{-10,180},{10,200}}),
              iconTransformation(
                extent={{-10,180},{10,200}})));

        BondGraph.BG_linear.Sensors.Sensor_p_e sensor_p_e1
          annotation (
            Placement(
              transformation(
                extent={{-160,-10},{-140,10}})));

        BondGraph.BG_linear.Sensors.Sensor_p_e sensor_p_e2
          annotation (
            Placement(
              transformation(
                extent={{-160,-50},{-140,-30}})));

        BondGraph.BG_linear.Sensors.Sensor_p_e sensor_p_e3
          annotation (
            Placement(
              transformation(
                extent={{140,-10},{160,10}})));

        BondGraph.BG_linear.Elements.Jzero jzero1(p=1, n=1)
          annotation (
            Placement(
              transformation(
                extent={{-120,-10},{-100,10}})));

        BondGraph.BG_linear.Elements.Jzero jzero2(p=1, n=1)
          annotation (
            Placement(
              transformation(
                extent={{-120,-50},{-100,-30}})));

        BondGraph.BG_linear.Elements.Jzero jzero3(n=1, p=2)
          annotation (
            Placement(
              transformation(
                extent={{100,-10},{120,10}})));

        BondGraph.Blocks.Const const1(c={{par_c_2}})
          annotation (
            Placement(
              transformation(
                extent={{-50,180},{-30,200}})));

        BondGraph.Blocks.Prop_sat prop_sat1(
          b={cat(   1,
                    {1},
                    par_c_1)},
          x_min={par_stat_1[1]},
          x_max={par_stat_1[size(par_stat_1, 1)]})
          annotation (
            Placement(
              transformation(
                extent={{-10,120},{10,140}})));

        BondGraph.Blocks.Prop_sat prop_sat2(
          b=identity(size(par_stat_1, 1)),
          x_min=1e-9*ones(size(par_stat_1, 1)),
          x_max=ones(size(par_stat_1, 1)))
          annotation (
            Placement(
              transformation(
                extent={{-10,30},{10,50}})));

        BondGraph.Blocks.Lin_sat lin_sat1(
          b={{par_dyn_1}},
          a={{par_dyn_2}},
          x_0={par_dyn_3},
          x_min={par_stat_1[1]},
          x_max={par_stat_1[size(par_stat_1, 1)]},
          x_der_min={par_dyn_4},
          x_der_max={par_dyn_5})
          annotation (
            Placement(
              transformation(
                extent={{-10,90},{10,110}})));

        BondGraph.Blocks.Sig_c sig_c1(
          p_x=par_stat_1,
          p_cn=par_stat_2,
          p_cp=par_stat_3)
          annotation (
            Placement(
              transformation(
                extent={{-10,60},{10,80}})));

        BondGraph.BG_nonlinear.Hydraulics.MHRT mhrt1(
          temp=par_res_1[1],
          part_mass_air=par_res_2[1],
          e_nom=par_res_3[1],
          f_nom=par_res_4[1],
          e_ref=par_res_5[1],
          f_ref=par_res_6[1],
          rho_ref=par_res_7[1],
          p=par_res_8[1],
          r_min=par_res_9[1],
          r_max=par_res_10[1],
          par_caus=par_res_11[1])
          annotation (
            Placement(
              transformation(
                extent={{-10,-10},{10,10}})));

        BondGraph.BG_nonlinear.Hydraulics.MHRT mhrt2(
          temp=par_res_1[2],
          part_mass_air=par_res_2[2],
          e_nom=par_res_3[2],
          f_nom=par_res_4[2],
          e_ref=par_res_5[2],
          f_ref=par_res_6[2],
          rho_ref=par_res_7[2],
          p=par_res_8[2],
          r_min=par_res_9[2],
          r_max=par_res_10[2],
          par_caus=par_res_11[2])
          annotation (
            Placement(
              transformation(
                extent={{-10,-50},{10,-30}})));

      equation
        connect(const1.port_out[1, 1], prop_sat1.port_in[1]) annotation (Line(
            points={{-40,181},{-40,170},{0,170},{0,139}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(port_in, prop_sat1.port_in[2]) annotation (Line(
            points={{0,190},{0,139}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(prop_sat1.port_out, lin_sat1.port_in) annotation (Line(
            points={{0,121},{0,109}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(lin_sat1.port_out[1], sig_c1.port_in) annotation (Line(
            points={{0,91},{0,79}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sig_c1.port_out, prop_sat2.port_in) annotation (Line(
            points={{0,61},{0,49}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(port_p1, sensor_p_e1.port_p) annotation (Line(
            points={{-190,0},{-159,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_p_e3.port_n, port_n) annotation (Line(
            points={{159,0},{190,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_p_e1.port_n, jzero1.port_p[1]) annotation (Line(
            points={{-141,0},{-119,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero3.port_n[1], sensor_p_e3.port_p) annotation (Line(
            points={{119,0},{141,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_p_e1.port_out, prop_sat1.port_in[3]) annotation (Line(
            points={{-150,-9},{-150,-20},{-130,-20},{-130,170},{0,170},{0,139}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(port_p2, sensor_p_e2.port_p) annotation (Line(
            points={{-190,-40},{-159,-40}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_p_e2.port_n,jzero2. port_p[1]) annotation (Line(
            points={{-141,-40},{-119,-40}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_p_e2.port_out, prop_sat1.port_in[4]) annotation (Line(
            points={{-150,-49},{-150,-60},{-130,-60},{-130,170},{0,170},{0,139}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_p_e3.port_out, prop_sat1.port_in[5]) annotation (Line(
            points={{150,-9},{150,-20},{130,-20},{130,170},{0,170},{0,139}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero1.port_n[1], mhrt1.port_p) annotation (Line(
            points={{-101,0},{-9,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(mhrt1.port_n, jzero3.port_p[1]) annotation (Line(
            points={{9,0},{90,0},{90,-0.5},{101,-0.5}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(prop_sat2.port_out[1], mhrt1.port_in) annotation (Line(
            points={{0,31},{0,9}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(prop_sat2.port_out[2], mhrt2.port_in) annotation (Line(
            points={{0,31},{0,20},{20,20},{20,-20},{0,-20},{0,-31}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero2.port_n[1], mhrt2.port_p) annotation (Line(
            points={{-101,-40},{-9,-40}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(mhrt2.port_n, jzero3.port_p[2]) annotation (Line(
            points={{9,-40},{90,-40},{90,0.5},{101,0.5}},
            color={0,0,0},
            smooth=Smooth.None));

        annotation (
          Diagram(
            coordinateSystem(
              preserveAspectRatio=true,
              extent={{-200,-200},{200,200}}),
            graphics),
          Icon(
            coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-200},{100,200}}),
            graphics={
                Rectangle(
                extent={{-100,200},{100,-200}},
                pattern=LinePattern.None,
                fillColor={0,0,0},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}),
                Rectangle(
                extent={{-80,180},{80,-180}},
                pattern=LinePattern.None,
                fillColor={181,117,52},
                fillPattern=FillPattern.Solid),
                        Text(
                  extent={{0,30},{0,-30}},
                  lineColor={0,0,0},
                  lineThickness=1,
                  fillColor={51,249,65},
                  fillPattern=FillPattern.Solid,
                textString="switch
3x2")}));
      end Switch_3x2;

      model Switch_4x3

        parameter Real par_c_1[5] = (1/0.1e5)*{1,0,0,0,0}
          "control gains, {u, p_1, p_2, p_3, p_4}";
        parameter Real par_c_2 = 0
          "control signal for par_c_1*{u, p_1, p_2, p_3, p_4} = 0";
        parameter Real par_dyn_1 = 1e6 "control signal gain";
        parameter Real par_dyn_2 = -1e6 "state gain";
        parameter Real par_dyn_3 = 0 "initial state";
        parameter Real par_dyn_4 = -1/10e-3 "minimal der(state) value";
        parameter Real par_dyn_5 = 1/10e-3 "maximal der(state) value";
        parameter Real par_stat_1[3] = {-1,0,1} "positions assignment";
        parameter Real par_stat_2[3] = {0.5,0.5,0.5}
          "coverage parameter for state < p_x[i]";
        parameter Real par_stat_3[3] = {0.5,0.5,0.5}
          "coverage parameter for state > p_x[i]";
        parameter Real par_res_1[4] = 313.15*{1,1,1,1} "temperature";
        parameter Real par_res_2[4] = 1.393e-6*{1,1,1,1}
          "mass proportion of undissolved air to total mass of mixture";
        parameter Real par_res_3[4] = 1e6*{1,1,1,1} "nominal effort";
        parameter Real par_res_4[4] = 1e-3*{1,1,1,1} "nominal flow";
        parameter Real par_res_5[4] = 0.5e5*{1,1,1,1}
          "reference pressure difference";
        parameter Real par_res_6[4] = 0.5e-3*{1,1,1,1} "reference volume flow";
        parameter Real par_res_7[4] = 854.307*{1,1,1,1} "reference density";
        parameter Real par_res_8[4] = 2*{1,1,1,1}
          "volume flow exponent for resistance calculation";
        parameter Real par_res_9[4] = 1e-2*par_res_3 ./ par_res_4
          "minimal resistance";
        parameter Real par_res_10[4] = 1e6*par_res_3 ./ par_res_4
          "maximal resistance";
        parameter Integer par_res_11[4] = 2*{1,1,1,1}
          "constitutive equation causality, par_res_5=1(effort out), par_res_5=2(flow out)";

        BondGraph.Interfaces.Port_p port_p1
          annotation (
            Placement(
              transformation(
                extent={{-200,-10},{-180,10}}),
              iconTransformation(
                extent={{-100,90},{-80,110}})));

        BondGraph.Interfaces.Port_p port_p2
          annotation (
            Placement(
              transformation(
                extent={{-200,-50},{-180,-30}}),
              iconTransformation(
                extent={{-100,-110},{-80,-90}})));

        BondGraph.Interfaces.Port_n port_n1
          annotation (
            Placement(
              transformation(
                extent={{180,-10},{200,10}}),
              iconTransformation(
                extent={{80,90},{100,110}})));

        BondGraph.Interfaces.Port_n port_n2
          annotation (
            Placement(
              transformation(
                extent={{180,-50},{200,-30}}),
              iconTransformation(
                extent={{80,-110},{100,-90}})));

        BondGraph.Interfaces.Port_in port_in
          annotation (
            Placement(
              transformation(
                extent={{-10,180},{10,200}}),
              iconTransformation(
                extent={{-10,180},{10,200}})));

        BondGraph.BG_linear.Sensors.Sensor_p_e sensor_p_e1
          annotation (
            Placement(
              transformation(
                extent={{-160,-10},{-140,10}})));

        BondGraph.BG_linear.Sensors.Sensor_p_e sensor_p_e2
          annotation (
            Placement(
              transformation(
                extent={{-160,-50},{-140,-30}})));

        BondGraph.BG_linear.Sensors.Sensor_p_e sensor_p_e3
          annotation (
            Placement(
              transformation(
                extent={{140,-10},{160,10}})));

        BondGraph.BG_linear.Sensors.Sensor_p_e sensor_p_e4
          annotation (
            Placement(
              transformation(
                extent={{140,-50},{160,-30}})));

        BondGraph.BG_linear.Elements.Jzero jzero1(p=1, n=2)
          annotation (
            Placement(
              transformation(
                extent={{-120,-10},{-100,10}})));

        BondGraph.BG_linear.Elements.Jzero jzero2(p=1, n=2)
          annotation (
            Placement(
              transformation(
                extent={{-120,-50},{-100,-30}})));

        BondGraph.BG_linear.Elements.Jzero jzero3(n=1, p=2)
          annotation (
            Placement(
              transformation(
                extent={{100,-10},{120,10}})));

        BondGraph.BG_linear.Elements.Jzero jzero4(n=1, p=2)
          annotation (
            Placement(
              transformation(
                extent={{100,-50},{120,-30}})));

        BondGraph.Blocks.Const const1(c={{par_c_2}})
          annotation (
            Placement(
              transformation(
                extent={{-50,180},{-30,200}})));

        BondGraph.Blocks.Prop_sat prop_sat1(
          b={cat(   1,
                    {1},
                    par_c_1)},
          x_min={par_stat_1[1]},
          x_max={par_stat_1[size(par_stat_1, 1)]})
          annotation (
            Placement(
              transformation(
                extent={{-10,120},{10,140}})));

        BondGraph.Blocks.Prop_sat prop_sat2(
          b=identity(size(par_stat_1, 1)),
          x_min=1e-9*ones(size(par_stat_1, 1)),
          x_max=ones(size(par_stat_1, 1)))
          annotation (
            Placement(
              transformation(
                extent={{-10,30},{10,50}})));

        BondGraph.Blocks.Lin_sat lin_sat1(
          b={{par_dyn_1}},
          a={{par_dyn_2}},
          x_0={par_dyn_3},
          x_min={par_stat_1[1]},
          x_max={par_stat_1[size(par_stat_1, 1)]},
          x_der_min={par_dyn_4},
          x_der_max={par_dyn_5})
          annotation (
            Placement(
              transformation(
                extent={{-10,90},{10,110}})));

        BondGraph.Blocks.Sig_c sig_c1(
          p_x=par_stat_1,
          p_cn=par_stat_2,
          p_cp=par_stat_3)
          annotation (
            Placement(
              transformation(
                extent={{-10,60},{10,80}})));

        BondGraph.BG_nonlinear.Hydraulics.MHRT mhrt1(
          temp=par_res_1[1],
          part_mass_air=par_res_2[1],
          e_nom=par_res_3[1],
          f_nom=par_res_4[1],
          e_ref=par_res_5[1],
          f_ref=par_res_6[1],
          rho_ref=par_res_7[1],
          p=par_res_8[1],
          r_min=par_res_9[1],
          r_max=par_res_10[1],
          par_caus=par_res_11[1])
          annotation (
            Placement(
              transformation(
                extent={{-10,-10},{10,10}})));

        BondGraph.BG_nonlinear.Hydraulics.MHRT mhrt2(
          temp=par_res_1[2],
          part_mass_air=par_res_2[2],
          e_nom=par_res_3[2],
          f_nom=par_res_4[2],
          e_ref=par_res_5[2],
          f_ref=par_res_6[2],
          rho_ref=par_res_7[2],
          p=par_res_8[2],
          r_min=par_res_9[2],
          r_max=par_res_10[2],
          par_caus=par_res_11[2])
          annotation (
            Placement(
              transformation(
                extent={{-10,-50},{10,-30}})));

        BondGraph.BG_nonlinear.Hydraulics.MHRT mhrt3(
          temp=par_res_1[3],
          part_mass_air=par_res_2[3],
          e_nom=par_res_3[3],
          f_nom=par_res_4[3],
          e_ref=par_res_5[3],
          f_ref=par_res_6[3],
          rho_ref=par_res_7[3],
          p=par_res_8[3],
          r_min=par_res_9[3],
          r_max=par_res_10[3],
          par_caus=par_res_11[3])
          annotation (
            Placement(
              transformation(
                extent={{-10,-90},{10,-70}})));

        BondGraph.BG_nonlinear.Hydraulics.MHRT mhrt4(
          temp=par_res_1[4],
          part_mass_air=par_res_2[4],
          e_nom=par_res_3[4],
          f_nom=par_res_4[4],
          e_ref=par_res_5[4],
          f_ref=par_res_6[4],
          rho_ref=par_res_7[4],
          p=par_res_8[4],
          r_min=par_res_9[4],
          r_max=par_res_10[4],
          par_caus=par_res_11[4])
          annotation (
            Placement(
              transformation(
                extent={{-10,-130},{10,-110}})));

      equation
        connect(const1.port_out[1, 1], prop_sat1.port_in[1]) annotation (Line(
            points={{-40,181},{-40,170},{0,170},{0,139}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(port_in, prop_sat1.port_in[2]) annotation (Line(
            points={{0,190},{0,139}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(prop_sat1.port_out, lin_sat1.port_in) annotation (Line(
            points={{0,121},{0,109}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(lin_sat1.port_out[1], sig_c1.port_in) annotation (Line(
            points={{0,91},{0,79}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sig_c1.port_out, prop_sat2.port_in) annotation (Line(
            points={{0,61},{0,49}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(port_p1, sensor_p_e1.port_p) annotation (Line(
            points={{-190,0},{-159,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_p_e3.port_n, port_n1) annotation (Line(
            points={{159,0},{190,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_p_e1.port_n, jzero1.port_p[1]) annotation (Line(
            points={{-141,0},{-119,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero3.port_n[1], sensor_p_e3.port_p) annotation (Line(
            points={{119,0},{141,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_p_e1.port_out, prop_sat1.port_in[3]) annotation (Line(
            points={{-150,-9},{-150,-20},{-130,-20},{-130,170},{0,170},{0,139}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(port_p2, sensor_p_e2.port_p) annotation (Line(
            points={{-190,-40},{-159,-40}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_p_e2.port_n,jzero2. port_p[1]) annotation (Line(
            points={{-141,-40},{-119,-40}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_p_e2.port_out, prop_sat1.port_in[4]) annotation (Line(
            points={{-150,-49},{-150,-60},{-130,-60},{-130,170},{0,170},{0,139}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_p_e3.port_out, prop_sat1.port_in[5]) annotation (Line(
            points={{150,-9},{150,-20},{130,-20},{130,170},{0,170},{0,139}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero4.port_n[1],sensor_p_e4. port_p) annotation (Line(
            points={{119,-40},{141,-40}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_p_e4.port_n,port_n2) annotation (Line(
            points={{159,-40},{190,-40}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_p_e4.port_out, prop_sat1.port_in[6]) annotation (Line(
            points={{150,-49},{150,-60},{130,-60},{130,170},{0,170},{0,139}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero1.port_n[1], mhrt1.port_p) annotation (Line(
            points={{-101,-0.5},{-80,-0.5},{-80,0},{-9,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero2.port_n[1], mhrt2.port_p) annotation (Line(
            points={{-101,-40.5},{-90,-40.5},{-90,-40},{-9,-40}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero1.port_n[2], mhrt3.port_p) annotation (Line(
            points={{-101,0.5},{-80,0.5},{-80,-80},{-9,-80}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero2.port_n[2], mhrt4.port_p) annotation (Line(
            points={{-101,-39.5},{-90,-39.5},{-90,-120},{-9,-120}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(mhrt1.port_n, jzero4.port_p[1]) annotation (Line(
            points={{9,0},{90,0},{90,-40.5},{101,-40.5}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(mhrt2.port_n, jzero3.port_p[1]) annotation (Line(
            points={{9,-40},{80,-40},{80,-0.5},{101,-0.5}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(mhrt3.port_n, jzero3.port_p[2]) annotation (Line(
            points={{9,-80},{80,-80},{80,0.5},{101,0.5}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(mhrt4.port_n, jzero4.port_p[2]) annotation (Line(
            points={{9,-120},{90,-120},{90,-39.5},{101,-39.5}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(prop_sat2.port_out[1], mhrt1.port_in) annotation (Line(
            points={{0,31},{0,9}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(prop_sat2.port_out[1], mhrt2.port_in) annotation (Line(
            points={{0,31},{0,20},{20,20},{20,-20},{0,-20},{0,-31}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(prop_sat2.port_out[3], mhrt3.port_in) annotation (Line(
            points={{0,31},{0,20},{20,20},{20,-60},{0,-60},{0,-71}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(prop_sat2.port_out[3], mhrt4.port_in) annotation (Line(
            points={{0,31},{0,20},{20,20},{20,-100},{0,-100},{0,-111}},
            color={0,0,0},
            smooth=Smooth.None));

        annotation (
          Diagram(
            coordinateSystem(
              preserveAspectRatio=true,
              extent={{-200,-200},{200,200}}),
            graphics),
          Icon(
            coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-200},{100,200}}),
            graphics={
                Rectangle(
                  extent={{-100,200},{100,-200}},
                  pattern=LinePattern.None,
                  fillColor={0,0,0},
                  fillPattern=FillPattern.Solid,
                lineColor={0,0,0}),
                Rectangle(
                extent={{-80,180},{80,-180}},
                pattern=LinePattern.None,
                fillColor={181,117,52},
                fillPattern=FillPattern.Solid),
                        Text(
                  extent={{0,30},{0,-30}},
                  lineColor={0,0,0},
                  lineThickness=1,
                  fillColor={51,249,65},
                  fillPattern=FillPattern.Solid,
                textString="switch
4x3")}));
      end Switch_4x3;
    end Valves;

    package Pipes
      model Pipe_x2

          parameter Real temp_p = 313.15 "temperature";
          parameter Real part_mass_air_p = 1.393e-6
          "mass proportion of undissolved air to total mass of mixture";
          parameter Real e_nom_p = 1e6 "nominal effort";
          parameter Real f_nom_p = 1e-3 "nominal flow";
          parameter Real r_min_p = 1e-2*e_nom_p/f_nom_p "minimal resistance";
          parameter Real r_max_p = 1e6*e_nom_p/f_nom_p "maximal resistance";
          parameter Real l_p = 1 "length";
          parameter Real d_p = 0.01 "hydraulic diameter";
          parameter Real c_p = 1e-8 "compressibility";
          parameter Real r_p = 1e-6 "relative hydraulic roughness";
          parameter Real e_0_p = 1e5 "initial pressure";
          parameter Real f_0_p = 0 "initial flow";
          parameter Integer hr_par_caus[2] = {2,2}
          "constitutive equation causality of hr, par_caus=1(effort out), par_caus=2(flow out)";
          parameter Real hrt_e_ref[2] = 0.5e5*{1,1}
          "reference pressure difference";
          parameter Real hrt_f_ref[2] = 0.5e-3*{1,1} "reference volume flow";
          parameter Real hrt_rho_ref[2] = 854.307*{1,1} "reference density";
          parameter Real hrt_p[2] = 2*{1,1}
          "volume flow exponent for resistance calculation";
          parameter Integer hrt_par_caus[2] = {2,2}
          "constitutive equation causality of hrt, par_caus=1(effort out), par_caus=2(flow out)";

          BondGraph.Interfaces.Port_p port_p
            annotation (
              Placement(
                transformation(
                  extent={{-300,-10},{-280,10}}),
                iconTransformation(
                  extent={{-300,-10},{-280,10}})));
          BondGraph.Interfaces.Port_n port_n
            annotation (
              Placement(
                transformation(
                  extent={{280,-10},{300,10}}),
                iconTransformation(
                  extent={{280,-10},{300,10}})));
          BondGraph.BG_linear.Elements.Jzero jzero1(p=1, n=1)
            annotation (
              Placement(
                transformation(
                  extent={{-260,-10},{-240,10}})));
          BondGraph.BG_linear.Elements.Jzero jzero2(p=1, n=1)
            annotation (
              Placement(
                transformation(
                  extent={{240,-10},{260,10}})));
          BondGraph.BG_nonlinear.Hydraulics.HRT hrt1(
          temp=temp_p,
          part_mass_air=part_mass_air_p,
          e_nom=e_nom_p,
          f_nom=f_nom_p,
          r_min=r_min_p,
          r_max=r_max_p,
          e_ref=hrt_e_ref[1],
          f_ref=hrt_f_ref[1],
          rho_ref=hrt_rho_ref[1],
          p=hrt_p[1],
          par_caus=hrt_par_caus[1])
            annotation (
              Placement(
                transformation(
                  extent={{-220,-10},{-200,10}})));
          BondGraph.BG_nonlinear.Hydraulics.HRT hrt2(
          temp=temp_p,
          part_mass_air=part_mass_air_p,
          e_nom=e_nom_p,
          f_nom=f_nom_p,
          r_min=r_min_p,
          r_max=r_max_p,
          e_ref=hrt_e_ref[2],
          f_ref=hrt_f_ref[2],
          rho_ref=hrt_rho_ref[2],
          p=hrt_p[2],
          par_caus=hrt_par_caus[2])
            annotation (
              Placement(
                transformation(
                  extent={{180,-10},{200,10}})));
          BondGraph.BG_nonlinear.Hydraulics.HR hr1(
          temp=temp_p,
          part_mass_air=part_mass_air_p,
          e_nom=e_nom_p,
          f_nom=f_nom_p,
          r_min=r_min_p,
          r_max=r_max_p,
          a=(0.5*asin(1))*d_p^2,
          l=(1/n_r)*l_p,
          d_h=d_p,
          r_h=r_p,
          par_caus=hr_par_caus[1])
            annotation (
              Placement(
                transformation(
                  extent={{-60,-10},{-40,10}})));
          BondGraph.BG_nonlinear.Hydraulics.HR hr2(
          temp=temp_p,
          part_mass_air=part_mass_air_p,
          e_nom=e_nom_p,
          f_nom=f_nom_p,
          r_min=r_min_p,
          r_max=r_max_p,
          a=(0.5*asin(1))*d_p^2,
          l=(1/n_r)*l_p,
          d_h=d_p,
          r_h=r_p,
          par_caus=hr_par_caus[2])
            annotation (
              Placement(
                transformation(
                  extent={{100,-10},{120,10}})));
          BondGraph.BG_nonlinear.Hydraulics.HC hc1(
          c=c_p,
          v_0=((1/n_c)*l_p)*(0.5*asin(1))*d_p^2,
          e_0=e_0_p,
          temp=temp_p,
          part_mass_air=part_mass_air_p,
          e_nom=e_nom_p,
          f_nom=f_nom_p)
            annotation (
              Placement(
                transformation(
                  extent={{-180,-10},{-160,10}})));
          BondGraph.BG_nonlinear.Hydraulics.HC hc2(
          c=c_p,
          v_0=((1/n_c)*l_p)*(0.5*asin(1))*d_p^2,
          e_0=e_0_p,
          temp=temp_p,
          part_mass_air=part_mass_air_p,
          e_nom=e_nom_p,
          f_nom=f_nom_p)
            annotation (
              Placement(
                transformation(
                  extent={{-100,-10},{-80,10}})));
          BondGraph.BG_nonlinear.Hydraulics.HC hc3(
          c=c_p,
          v_0=((1/n_c)*l_p)*(0.5*asin(1))*d_p^2,
          e_0=e_0_p,
          temp=temp_p,
          part_mass_air=part_mass_air_p,
          e_nom=e_nom_p,
          f_nom=f_nom_p)
            annotation (
              Placement(
                transformation(
                  extent={{-20,-10},{0,10}})));
          BondGraph.BG_nonlinear.Hydraulics.HC hc4(
          c=c_p,
          v_0=((1/n_c)*l_p)*(0.5*asin(1))*d_p^2,
          e_0=e_0_p,
          temp=temp_p,
          part_mass_air=part_mass_air_p,
          e_nom=e_nom_p,
          f_nom=f_nom_p)
            annotation (
              Placement(
                transformation(
                  extent={{60,-10},{80,10}})));
          BondGraph.BG_nonlinear.Hydraulics.HC hc5(
          c=c_p,
          v_0=((1/n_c)*l_p)*(0.5*asin(1))*d_p^2,
          e_0=e_0_p,
          temp=temp_p,
          part_mass_air=part_mass_air_p,
          e_nom=e_nom_p,
          f_nom=f_nom_p)
            annotation (
              Placement(
                transformation(
                  extent={{140,-10},{160,10}})));
          BondGraph.BG_nonlinear.Hydraulics.HI hi1(
          temp=temp_p,
          part_mass_air=part_mass_air_p,
          e_nom=e_nom_p,
          f_nom=f_nom_p,
          a=(0.5*asin(1))*d_p^2,
          l=(1/n_i)*l_p,
          f_0=f_0_p)
            annotation (
              Placement(
                transformation(
                  extent={{-140,-10},{-120,10}})));
          BondGraph.BG_nonlinear.Hydraulics.HI hi2(
          temp=temp_p,
          part_mass_air=part_mass_air_p,
          e_nom=e_nom_p,
          f_nom=f_nom_p,
          a=(0.5*asin(1))*d_p^2,
          l=(1/n_i)*l_p,
          f_0=f_0_p)
            annotation (
              Placement(
                transformation(
                  extent={{20,-10},{40,10}})));

      protected
          parameter Real n_i = 2;
          parameter Real n_r = n_i;
          parameter Real n_c = 5;

      equation
          connect(port_p, jzero1.port_p[1]) annotation (Line(
              points={{-290,0},{-259,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jzero2.port_n[1], port_n) annotation (Line(
              points={{259,0},{290,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jzero1.port_n[1], hrt1.port_p) annotation (Line(
              points={{-241,0},{-219,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(hrt1.port_n, hc1.port_p) annotation (Line(
              points={{-201,0},{-179,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(hc1.port_n, hi1.port_p) annotation (Line(
              points={{-161,0},{-139,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(hi1.port_n, hc2.port_p) annotation (Line(
              points={{-121,0},{-99,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(hc2.port_n, hr1.port_p) annotation (Line(
              points={{-81,0},{-59,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(hr1.port_n, hc3.port_p) annotation (Line(
              points={{-41,0},{-19,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(hc3.port_n, hi2.port_p) annotation (Line(
              points={{-1,0},{21,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(hi2.port_n, hc4.port_p) annotation (Line(
              points={{39,0},{61,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(hc4.port_n, hr2.port_p) annotation (Line(
              points={{79,0},{101,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(hr2.port_n, hc5.port_p) annotation (Line(
              points={{119,0},{141,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(hc5.port_n, hrt2.port_p) annotation (Line(
              points={{159,0},{181,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(hrt2.port_n, jzero2.port_p[1]) annotation (Line(
              points={{199,0},{241,0}},
              color={0,0,0},
              smooth=Smooth.None));

          annotation (
            Diagram(
              coordinateSystem(
                preserveAspectRatio=true,
                extent={{-300,-100},{300,100}}),
              graphics),
            Icon(
              coordinateSystem(
                preserveAspectRatio=true,
                extent={{-300,-100},{300,100}}),
              graphics={Rectangle(
                  extent={{-300,100},{300,80}},
                  fillColor={0,0,0},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None,
                  lineColor={0,0,0}),
                Rectangle(
                  extent={{-300,-80},{300,-100}},
                  fillColor={0,0,0},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None,
                  lineColor={0,0,0}),
                Rectangle(
                  extent={{-300,80},{300,-80}},
                  pattern=LinePattern.None,
                  fillColor={181,117,52},
                  fillPattern=FillPattern.Solid),
                        Text(
                  extent={{0,30},{0,-30}},
                  lineColor={0,0,0},
                  lineThickness=1,
                  fillColor={51,249,65},
                  fillPattern=FillPattern.Solid,
                textString="pipe_x2")}));
      end Pipe_x2;

        model Pipe_x2_lin

          parameter Real r_p[4] = 1e9*ones(4) "resistance parameters";
          parameter Real c_p[5] = 1e-12*ones(5) "compliance parameters";
          parameter Real i_p[2] = 1e-1*ones(2) "inductance parameters";
          parameter Real e_0_p = 1e5 "initial pressure";
          parameter Real f_0_p = 0 "initial flow";

          BondGraph.Interfaces.Port_p port_p
            annotation (
              Placement(
                transformation(
                  extent={{-300,-10},{-280,10}}),
                iconTransformation(
                  extent={{-300,-10},{-280,10}})));
          BondGraph.Interfaces.Port_n port_n
            annotation (
              Placement(
                transformation(
                  extent={{280,-10},{300,10}}),
                iconTransformation(
                  extent={{280,-10},{300,10}})));
          BondGraph.BG_linear.Elements.Jone jone1(p=1, n=2)
            annotation (
              Placement(
                transformation(
                  extent={{-220,-10},{-200,10}})));
          BondGraph.BG_linear.Elements.Jone jone2(p=1, n=2)
            annotation (
              Placement(
                transformation(
                  extent={{-140,-10},{-120,10}})));
          BondGraph.BG_linear.Elements.Jone jone3(p=1, n=2)
            annotation (
              Placement(
                transformation(
                  extent={{-60,-10},{-40,10}})));
          BondGraph.BG_linear.Elements.Jone jone4(p=1, n=2)
            annotation (
              Placement(
                transformation(
                  extent={{20,-10},{40,10}})));
          BondGraph.BG_linear.Elements.Jone jone5(p=1, n=2)
            annotation (
              Placement(
                transformation(
                  extent={{100,-10},{120,10}})));
          BondGraph.BG_linear.Elements.Jone jone6(p=1, n=2)
            annotation (
              Placement(
                transformation(
                  extent={{180,-10},{200,10}})));
          BondGraph.BG_linear.Elements.Jzero jzero1(p=1, n=1)
            annotation (
              Placement(
                transformation(
                  extent={{-260,-10},{-240,10}})));
          BondGraph.BG_linear.Elements.Jzero jzero2(p=1, n=1)
            annotation (
              Placement(
                transformation(
                  extent={{240,-10},{260,10}})));
          BondGraph.BG_linear.Elements.Jzero jzero3(p=1, n=2)
            annotation (
              Placement(
                transformation(
                  extent={{-180,-10},{-160,10}})));
          BondGraph.BG_linear.Elements.Jzero jzero4(p=1, n=2)
            annotation (
              Placement(
                transformation(
                  extent={{-100,-10},{-80,10}})));
          BondGraph.BG_linear.Elements.Jzero jzero5(p=1, n=2)
            annotation (
              Placement(
                transformation(
                  extent={{-20,-10},{0,10}})));
          BondGraph.BG_linear.Elements.Jzero jzero6(p=1, n=2)
            annotation (
              Placement(
                transformation(
                  extent={{60,-10},{80,10}})));
          BondGraph.BG_linear.Elements.Jzero jzero7(p=1, n=2)
            annotation (
              Placement(
                transformation(
                  extent={{140,-10},{160,10}})));
          BondGraph.BG_linear.Elements.R r1(r=r_p[1])
            annotation (
              Placement(
                transformation(
                  extent={{-220,30},{-200,50}})));
          BondGraph.BG_linear.Elements.R r2(r=r_p[2])
            annotation (
              Placement(
                transformation(
                  extent={{-60,30},{-40,50}})));
          BondGraph.BG_linear.Elements.R r3(r=r_p[3])
            annotation (
              Placement(
                transformation(
                  extent={{100,30},{120,50}})));
          BondGraph.BG_linear.Elements.R r4(r=r_p[4])
            annotation (
              Placement(
                transformation(
                  extent={{180,30},{200,50}})));
          BondGraph.BG_linear.Elements.C c1(c=c_p[1], e_0=e_0_p)
            annotation (
              Placement(
                transformation(
                  extent={{-180,30},{-160,50}})));
          BondGraph.BG_linear.Elements.C c2(c=c_p[2], e_0=e_0_p)
            annotation (
              Placement(
                transformation(
                  extent={{-100,30},{-80,50}})));
          BondGraph.BG_linear.Elements.C c3(c=c_p[3], e_0=e_0_p)
            annotation (
              Placement(
                transformation(
                  extent={{-20,30},{0,50}})));
          BondGraph.BG_linear.Elements.C c4(c=c_p[4], e_0=e_0_p)
            annotation (
              Placement(
                transformation(
                  extent={{60,30},{80,50}})));
          BondGraph.BG_linear.Elements.C c5(c=c_p[5], e_0=e_0_p)
            annotation (
              Placement(
                transformation(
                  extent={{140,30},{160,50}})));
          BondGraph.BG_linear.Elements.I i1(i=i_p[1], f_0=f_0_p)
            annotation (
              Placement(
                transformation(
                  extent={{-140,30},{-120,50}})));
          BondGraph.BG_linear.Elements.I i2(i=i_p[2], f_0=f_0_p)
            annotation (
              Placement(
                transformation(
                  extent={{20,30},{40,50}})));

        equation
          connect(port_p, jzero1.port_p[1]) annotation (Line(
              points={{-290,0},{-259,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jzero2.port_n[1], port_n) annotation (Line(
              points={{259,0},{290,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jone1.port_n[2], r1.port) annotation (Line(
              points={{-201,0.5},{-201,0.5},{-201,40}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jzero3.port_n[2], c1.port) annotation (Line(
              points={{-161,0.5},{-161,0.5},{-161,40}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jone2.port_n[2], i1.port) annotation (Line(
              points={{-121,0.5},{-121,0.5},{-121,40}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jzero4.port_n[2], c2.port) annotation (Line(
              points={{-81,0.5},{-81,40}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jone1.port_n[1], jzero3.port_p[1]) annotation (Line(
              points={{-201,-0.5},{-189.5,-0.5},{-189.5,0},{-179,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jzero3.port_n[1], jone2.port_p[1]) annotation (Line(
              points={{-161,-0.5},{-149.5,-0.5},{-149.5,0},{-139,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jone2.port_n[1], jzero4.port_p[1]) annotation (Line(
              points={{-121,-0.5},{-109.5,-0.5},{-109.5,0},{-99,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jone3.port_n[2],r2. port) annotation (Line(
              points={{-41,0.5},{-41,40}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jzero5.port_n[2],c3. port) annotation (Line(
              points={{-1,0.5},{-1,40}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jone4.port_n[2],i2. port) annotation (Line(
              points={{39,0.5},{39,40}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jzero6.port_n[2],c4. port) annotation (Line(
              points={{79,0.5},{79,40}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jone3.port_n[1],jzero5. port_p[1]) annotation (Line(
              points={{-41,-0.5},{-29.5,-0.5},{-29.5,0},{-19,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jzero5.port_n[1],jone4. port_p[1]) annotation (Line(
              points={{-1,-0.5},{10.5,-0.5},{10.5,0},{21,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jone4.port_n[1],jzero6. port_p[1]) annotation (Line(
              points={{39,-0.5},{50.5,-0.5},{50.5,0},{61,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jone5.port_n[2],r3. port) annotation (Line(
              points={{119,0.5},{119,40}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jone5.port_n[1],jzero7. port_p[1]) annotation (Line(
              points={{119,-0.5},{130.5,-0.5},{130.5,0},{141,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jzero7.port_n[2],c5. port) annotation (Line(
              points={{159,0.5},{159,40}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jone6.port_n[2],r4. port) annotation (Line(
              points={{199,0.5},{199,40}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jzero4.port_n[1], jone3.port_p[1]) annotation (Line(
              points={{-81,-0.5},{-69.5,-0.5},{-69.5,0},{-59,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jzero6.port_n[1], jone5.port_p[1]) annotation (Line(
              points={{79,-0.5},{90.5,-0.5},{90.5,0},{101,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jzero7.port_n[1], jone6.port_p[1]) annotation (Line(
              points={{159,-0.5},{170.5,-0.5},{170.5,0},{181,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jone6.port_n[1], jzero2.port_p[1]) annotation (Line(
              points={{199,-0.5},{220.5,-0.5},{220.5,0},{241,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jzero1.port_n[1], jone1.port_p[1]) annotation (Line(
              points={{-241,0},{-219,0}},
              color={0,0,0},
              smooth=Smooth.None));

          annotation (
            Diagram(
              coordinateSystem(
                preserveAspectRatio=true,
                extent={{-300,-100},{300,100}}),
              graphics),
            Icon(
              coordinateSystem(
                preserveAspectRatio=true,
                extent={{-300,-100},{300,100}}),
              graphics={Rectangle(
                  extent={{-300,100},{300,80}},
                  fillColor={0,0,0},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None,
                  lineColor={0,0,0}),
                Rectangle(
                  extent={{-300,-80},{300,-100}},
                  fillColor={0,0,0},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None,
                  lineColor={0,0,0}),
                Rectangle(
                  extent={{-300,80},{300,-80}},
                  pattern=LinePattern.None,
                  fillColor={125,190,255},
                  fillPattern=FillPattern.Solid),
                        Text(
                  extent={{0,30},{0,-30}},
                  lineColor={0,0,0},
                  lineThickness=1,
                  fillColor={51,249,65},
                  fillPattern=FillPattern.Solid,
                  textString="pipe_x2_lin")}));
        end Pipe_x2_lin;

      model Pipe_x1

          parameter Real temp_p = 313.15 "temperature";
          parameter Real part_mass_air_p = 1.393e-6
          "mass proportion of undissolved air to total mass of mixture";
          parameter Real e_nom_p = 1e6 "nominal effort";
          parameter Real f_nom_p = 1e-3 "nominal flow";
          parameter Real r_min_p = 1e-2*e_nom_p/f_nom_p "minimal resistance";
          parameter Real r_max_p = 1e6*e_nom_p/f_nom_p "maximal resistance";
          parameter Real l_p = 1 "length";
          parameter Real d_p = 0.01 "hydraulic diameter";
          parameter Real c_p = 1e-8 "compressibility";
          parameter Real r_p = 1e-6 "relative hydraulic roughness";
          parameter Real e_0_p = 1e5 "initial pressure";
          parameter Real f_0_p = 0 "initial flow";
          parameter Integer hr_par_caus[2] = {2,2}
          "constitutive equation causality of hr, par_caus=1(effort out), par_caus=2(flow out)";
          parameter Real hrt_e_ref[2] = 0.5e5*{1,1}
          "reference pressure difference";
          parameter Real hrt_f_ref[2] = 0.5e-3*{1,1} "reference volume flow";
          parameter Real hrt_rho_ref[2] = 854.307*{1,1} "reference density";
          parameter Real hrt_p[2] = 2*{1,1}
          "volume flow exponent for resistance calculation";
          parameter Integer hrt_par_caus[2] = {2,2}
          "constitutive equation causality of hrt, par_caus=1(effort out), par_caus=2(flow out)";

          BondGraph.Interfaces.Port_p port_p
            annotation (
              Placement(
                transformation(
                  extent={{-300,-10},{-280,10}}),
                iconTransformation(
                  extent={{-300,-10},{-280,10}})));
          BondGraph.Interfaces.Port_n port_n
            annotation (
              Placement(
                transformation(
                  extent={{280,-10},{300,10}}),
                iconTransformation(
                  extent={{280,-10},{300,10}})));
          BondGraph.BG_linear.Elements.Jzero jzero1(p=1, n=1)
            annotation (
              Placement(
                transformation(
                  extent={{-260,-10},{-240,10}})));
          BondGraph.BG_linear.Elements.Jzero jzero2(p=1, n=1)
            annotation (
              Placement(
                transformation(
                  extent={{240,-10},{260,10}})));
          BondGraph.BG_nonlinear.Hydraulics.HRT hrt1(
          temp=temp_p,
          part_mass_air=part_mass_air_p,
          e_nom=e_nom_p,
          f_nom=f_nom_p,
          r_min=r_min_p,
          r_max=r_max_p,
          e_ref=hrt_e_ref[1],
          f_ref=hrt_f_ref[1],
          rho_ref=hrt_rho_ref[1],
          p=hrt_p[1],
          par_caus=hrt_par_caus[1])
            annotation (
              Placement(
                transformation(
                  extent={{-220,-10},{-200,10}})));
          BondGraph.BG_nonlinear.Hydraulics.HRT hrt2(
          temp=temp_p,
          part_mass_air=part_mass_air_p,
          e_nom=e_nom_p,
          f_nom=f_nom_p,
          r_min=r_min_p,
          r_max=r_max_p,
          e_ref=hrt_e_ref[2],
          f_ref=hrt_f_ref[2],
          rho_ref=hrt_rho_ref[2],
          p=hrt_p[2],
          par_caus=hrt_par_caus[2])
            annotation (
              Placement(
                transformation(
                  extent={{20,-10},{40,10}})));
          BondGraph.BG_nonlinear.Hydraulics.HR hr1(
          temp=temp_p,
          part_mass_air=part_mass_air_p,
          e_nom=e_nom_p,
          f_nom=f_nom_p,
          r_min=r_min_p,
          r_max=r_max_p,
          a=(0.5*asin(1))*d_p^2,
          l=(1/n_r)*l_p,
          d_h=d_p,
          r_h=r_p,
          par_caus=hr_par_caus[1])
            annotation (
              Placement(
                transformation(
                  extent={{-60,-10},{-40,10}})));
          BondGraph.BG_nonlinear.Hydraulics.HC hc1(
          c=c_p,
          v_0=((1/n_c)*l_p)*(0.5*asin(1))*d_p^2,
          e_0=e_0_p,
          temp=temp_p,
          part_mass_air=part_mass_air_p,
          e_nom=e_nom_p,
          f_nom=f_nom_p)
            annotation (
              Placement(
                transformation(
                  extent={{-180,-10},{-160,10}})));
          BondGraph.BG_nonlinear.Hydraulics.HC hc2(
          c=c_p,
          v_0=((1/n_c)*l_p)*(0.5*asin(1))*d_p^2,
          e_0=e_0_p,
          temp=temp_p,
          part_mass_air=part_mass_air_p,
          e_nom=e_nom_p,
          f_nom=f_nom_p)
            annotation (
              Placement(
                transformation(
                  extent={{-100,-10},{-80,10}})));
          BondGraph.BG_nonlinear.Hydraulics.HC hc3(
          c=c_p,
          v_0=((1/n_c)*l_p)*(0.5*asin(1))*d_p^2,
          e_0=e_0_p,
          temp=temp_p,
          part_mass_air=part_mass_air_p,
          e_nom=e_nom_p,
          f_nom=f_nom_p)
            annotation (
              Placement(
                transformation(
                  extent={{-20,-10},{0,10}})));
          BondGraph.BG_nonlinear.Hydraulics.HI hi1(
          temp=temp_p,
          part_mass_air=part_mass_air_p,
          e_nom=e_nom_p,
          f_nom=f_nom_p,
          a=(0.5*asin(1))*d_p^2,
          l=(1/n_i)*l_p,
          f_0=f_0_p)
            annotation (
              Placement(
                transformation(
                  extent={{-140,-10},{-120,10}})));

      protected
          parameter Real n_i = 1;
          parameter Real n_r = n_i;
          parameter Real n_c = 3;

      equation
          connect(port_p, jzero1.port_p[1]) annotation (Line(
              points={{-290,0},{-259,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jzero2.port_n[1], port_n) annotation (Line(
              points={{259,0},{290,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jzero1.port_n[1], hrt1.port_p) annotation (Line(
              points={{-241,0},{-219,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(hrt1.port_n, hc1.port_p) annotation (Line(
              points={{-201,0},{-179,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(hc1.port_n, hi1.port_p) annotation (Line(
              points={{-161,0},{-139,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(hi1.port_n, hc2.port_p) annotation (Line(
              points={{-121,0},{-99,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(hc2.port_n, hr1.port_p) annotation (Line(
              points={{-81,0},{-59,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(hr1.port_n, hc3.port_p) annotation (Line(
              points={{-41,0},{-19,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(hrt2.port_n, jzero2.port_p[1]) annotation (Line(
              points={{39,0},{241,0}},
              color={0,0,0},
              smooth=Smooth.None));

        connect(hc3.port_n, hrt2.port_p) annotation (Line(
            points={{-1,0},{21,0}},
            color={0,0,0},
            smooth=Smooth.None));
          annotation (
            Diagram(
              coordinateSystem(
                preserveAspectRatio=true,
                extent={{-300,-100},{300,100}}),
              graphics),
            Icon(
              coordinateSystem(
                preserveAspectRatio=true,
                extent={{-300,-100},{300,100}}),
              graphics={Rectangle(
                  extent={{-300,100},{300,80}},
                  fillColor={0,0,0},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None,
                  lineColor={0,0,0}),
                Rectangle(
                  extent={{-300,-80},{300,-100}},
                  fillColor={0,0,0},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None,
                  lineColor={0,0,0}),
                Rectangle(
                  extent={{-300,80},{300,-80}},
                  pattern=LinePattern.None,
                  fillColor={181,117,52},
                  fillPattern=FillPattern.Solid),
                        Text(
                  extent={{0,30},{0,-30}},
                  lineColor={0,0,0},
                  lineThickness=1,
                  fillColor={51,249,65},
                  fillPattern=FillPattern.Solid,
                textString="pipe_x1")}));
      end Pipe_x1;

        model Pipe_x1_lin

          parameter Real r_p[3] = 1e9*ones(3) "resistance parameters";
          parameter Real c_p[3] = 1e-12*ones(3) "compliance parameters";
          parameter Real i_p[1] = 1e-1*ones(1) "inductance parameters";
          parameter Real e_0_p = 1e5 "initial pressure";
          parameter Real f_0_p = 0 "initial flow";

          BondGraph.Interfaces.Port_p port_p
            annotation (
              Placement(
                transformation(
                  extent={{-300,-10},{-280,10}}),
                iconTransformation(
                  extent={{-300,-10},{-280,10}})));
          BondGraph.Interfaces.Port_n port_n
            annotation (
              Placement(
                transformation(
                  extent={{280,-10},{300,10}}),
                iconTransformation(
                  extent={{280,-10},{300,10}})));
          BondGraph.BG_linear.Elements.Jone jone1(p=1, n=2)
            annotation (
              Placement(
                transformation(
                  extent={{-220,-10},{-200,10}})));
          BondGraph.BG_linear.Elements.Jone jone2(p=1, n=2)
            annotation (
              Placement(
                transformation(
                  extent={{-140,-10},{-120,10}})));
          BondGraph.BG_linear.Elements.Jone jone3(p=1, n=2)
            annotation (
              Placement(
                transformation(
                  extent={{-60,-10},{-40,10}})));
          BondGraph.BG_linear.Elements.Jone jone6(n=2, p=1)
            annotation (
              Placement(
                transformation(
                  extent={{20,-10},{40,10}})));
          BondGraph.BG_linear.Elements.Jzero jzero1(p=1, n=1)
            annotation (
              Placement(
                transformation(
                  extent={{-260,-10},{-240,10}})));
          BondGraph.BG_linear.Elements.Jzero jzero2(p=1, n=1)
            annotation (
              Placement(
                transformation(
                  extent={{240,-10},{260,10}})));
          BondGraph.BG_linear.Elements.Jzero jzero3(p=1, n=2)
            annotation (
              Placement(
                transformation(
                  extent={{-180,-10},{-160,10}})));
          BondGraph.BG_linear.Elements.Jzero jzero4(p=1, n=2)
            annotation (
              Placement(
                transformation(
                  extent={{-100,-10},{-80,10}})));
          BondGraph.BG_linear.Elements.Jzero jzero5(p=1, n=2)
            annotation (
              Placement(
                transformation(
                  extent={{-20,-10},{0,10}})));
          BondGraph.BG_linear.Elements.R r1(r=r_p[1])
            annotation (
              Placement(
                transformation(
                  extent={{-220,20},{-200,40}})));
          BondGraph.BG_linear.Elements.R r2(r=r_p[2])
            annotation (
              Placement(
                transformation(
                  extent={{-60,20},{-40,40}})));
          BondGraph.BG_linear.Elements.R r4(r=r_p[3])
            annotation (
              Placement(
                transformation(
                  extent={{20,20},{40,40}})));
          BondGraph.BG_linear.Elements.C c1(c=c_p[1], e_0=e_0_p)
            annotation (
              Placement(
                transformation(
                  extent={{-180,20},{-160,40}})));
          BondGraph.BG_linear.Elements.C c2(c=c_p[2], e_0=e_0_p)
            annotation (
              Placement(
                transformation(
                  extent={{-100,20},{-80,40}})));
          BondGraph.BG_linear.Elements.C c3(c=c_p[3], e_0=e_0_p)
            annotation (
              Placement(
                transformation(
                  extent={{-20,20},{0,40}})));
          BondGraph.BG_linear.Elements.I i1(i=i_p[1], f_0=f_0_p)
            annotation (
              Placement(
                transformation(
                  extent={{-140,20},{-120,40}})));

        equation
          connect(port_p, jzero1.port_p[1]) annotation (Line(
              points={{-290,0},{-259,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jzero2.port_n[1], port_n) annotation (Line(
              points={{259,0},{290,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jone1.port_n[2], r1.port) annotation (Line(
              points={{-201,0.5},{-201,30}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jzero3.port_n[2], c1.port) annotation (Line(
              points={{-161,0.5},{-161,30}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jone2.port_n[2], i1.port) annotation (Line(
              points={{-121,0.5},{-121,30}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jzero4.port_n[2], c2.port) annotation (Line(
              points={{-81,0.5},{-81,30}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jone1.port_n[1], jzero3.port_p[1]) annotation (Line(
              points={{-201,-0.5},{-189.5,-0.5},{-189.5,0},{-179,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jzero3.port_n[1], jone2.port_p[1]) annotation (Line(
              points={{-161,-0.5},{-149.5,-0.5},{-149.5,0},{-139,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jone2.port_n[1], jzero4.port_p[1]) annotation (Line(
              points={{-121,-0.5},{-109.5,-0.5},{-109.5,0},{-99,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jone3.port_n[2],r2. port) annotation (Line(
              points={{-41,0.5},{-41,30}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jone3.port_n[1],jzero5. port_p[1]) annotation (Line(
              points={{-41,-0.5},{-29.5,-0.5},{-29.5,0},{-19,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jone6.port_n[2],r4. port) annotation (Line(
              points={{39,0.5},{39,30}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jzero4.port_n[1], jone3.port_p[1]) annotation (Line(
              points={{-81,-0.5},{-69.5,-0.5},{-69.5,0},{-59,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jone6.port_n[1], jzero2.port_p[1]) annotation (Line(
              points={{39,-0.5},{60.5,-0.5},{60.5,0},{241,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(jzero1.port_n[1], jone1.port_p[1]) annotation (Line(
              points={{-241,0},{-219,0}},
              color={0,0,0},
              smooth=Smooth.None));

          connect(jzero5.port_n[1], jone6.port_p[1]) annotation (Line(
              points={{-1,-0.5},{9.5,-0.5},{9.5,0},{21,0}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(c3.port, jzero5.port_n[2]) annotation (Line(
              points={{-1,30},{-1,0.5}},
              color={0,0,0},
              smooth=Smooth.None));
          annotation (
            Diagram(
              coordinateSystem(
                preserveAspectRatio=true,
                extent={{-300,-100},{300,100}}),
              graphics),
            Icon(
              coordinateSystem(
                preserveAspectRatio=true,
                extent={{-300,-100},{300,100}}),
              graphics={Rectangle(
                  extent={{-300,100},{300,80}},
                  fillColor={0,0,0},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None,
                  lineColor={0,0,0}),
                Rectangle(
                  extent={{-300,-80},{300,-100}},
                  fillColor={0,0,0},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None,
                  lineColor={0,0,0}),
                Rectangle(
                  extent={{-300,80},{300,-80}},
                  pattern=LinePattern.None,
                  fillColor={125,190,255},
                  fillPattern=FillPattern.Solid),
                        Text(
                  extent={{0,30},{0,-30}},
                  lineColor={0,0,0},
                  lineThickness=1,
                  fillColor={51,249,65},
                  fillPattern=FillPattern.Solid,
                  textString="pipe_x1_lin")}));
        end Pipe_x1_lin;
    end Pipes;

    package Cylinders
      model Cylinder_lin
        parameter Real r_tf[2] = (1/(0.5*asin(1)*0.1^2))*{1,1}
          "transformer parameters";
        parameter Real r_hydr[2] = 1e8*{1,1} "hydraulic resistance parameters";
        parameter Real c_hydr[2] = 1e-12*{1,1}
          "hydraulic capacitance parameters";
        parameter Real e_0_hydr[2] = 1e5*{1,1} "initial pressures";
        parameter Real r_mech_1 = 1e3 "mechanical resistance parameter, piston";
        parameter Real i_mech = 1 "mechanical inductance parameter, piston";
        parameter Real f_0_mech = 0 "initial velocity, piston";
        parameter Real c_mech_lim = 1e-9
          "mechanical capacitance parameter, limiter";
        parameter Real q_mech_lim[2] = {0,1e-1} "position limits, piston";
        parameter Real q_0_mech = q_mech_lim[1] "initial position, piston";
        parameter Real c_mech = 1e-9
          "mechanical capacitance parameter, piston rod";
        parameter Real e_0_mech = 0 "initial force, piston rod";
        parameter Real r_mech_2 = 1e3
          "mechanical resistance parameter, piston rod";

        BondGraph.Interfaces.Port_p port_p1
          annotation (
            Placement(
              transformation(
                extent={{-200,70},{-180,90}}),
              iconTransformation(
                extent={{-100,90},{-80,110}})));

        BondGraph.Interfaces.Port_p port_p2
          annotation (
            Placement(
              transformation(
                extent={{-200,-90},{-180,-70}}),
              iconTransformation(
                extent={{-100,-110},{-80,-90}})));

        BondGraph.Interfaces.Port_n port_n
          annotation (
            Placement(
              transformation(
                extent={{180,-10},{200,10}}),
              iconTransformation(
                extent={{-10,-200},{10,-180}})));

        BondGraph.BG_linear.Elements.TF tf1(r=r_tf[1])
          annotation (
            Placement(
              transformation(
                extent={{-80,70},{-60,90}})));

        BondGraph.BG_linear.Elements.TF tf2(r=-r_tf[2])
          annotation (
            Placement(
              transformation(
                extent={{-80,-90},{-60,-70}})));

        BondGraph.BG_linear.Elements.C c_h_1(c=c_hydr[1], e_0=e_0_hydr[1])
          annotation (
            Placement(
              transformation(
                extent={{-120,40},{-100,60}})));

        BondGraph.BG_linear.Elements.C c_h_2(c=c_hydr[2], e_0=e_0_hydr[2])
          annotation (
            Placement(
              transformation(
                extent={{-120,-60},{-100,-40}})));

        BondGraph.BG_linear.Elements.R r_h_1(r=r_hydr[1])
          annotation (
            Placement(
              transformation(
                extent={{-160,40},{-140,60}})));

        BondGraph.BG_linear.Elements.R r_h_2(r=r_hydr[2])
          annotation (
            Placement(
              transformation(
                extent={{-160,-60},{-140,-40}})));

        BondGraph.BG_linear.Elements.I i_m_1(i=i_mech, f_0=f_0_mech)
          annotation (
            Placement(
              transformation(
                extent={{40,70},{60,90}})));

        BondGraph.BG_linear.Elements.R r_m_1(r=r_mech_1)
          annotation (
            Placement(
              transformation(
                extent={{40,40},{60,60}})));

        BondGraph.BG_linear.Elements.MSe mse1
          annotation (
            Placement(
              transformation(
                extent={{40,-20},{60,0}})));

        BondGraph.BG_linear.Sensors.Sensor_i_f sensor_i_f1(gain=1/c_mech_lim,
            y_0=q_0_mech)
          annotation (
            Placement(
              transformation(
                extent={{-40,-60},{-20,-40}})));

        Modelica.Blocks.Nonlinear.DeadZone deadZone(uMax=q_mech_lim[2]/c_mech_lim,
            uMin=q_mech_lim[1]/c_mech_lim)
          annotation (
            Placement(
              transformation(
                extent={{0,10},{20,30}})));

        BondGraph.BG_linear.Elements.C c_m_1(c=c_mech, e_0=e_0_mech)
          annotation (
            Placement(
              transformation(
                extent={{140,10},{160,30}})));

        BondGraph.BG_linear.Elements.R r_m_2(r=r_mech_2)
          annotation (
            Placement(
              transformation(
                extent={{140,70},{160,90}})));

        BondGraph.BG_linear.Elements.Jone jone1(p=1, n=2)
          annotation (
            Placement(
              transformation(
                  extent={{-160,70},{-140,90}})));

        BondGraph.BG_linear.Elements.Jone jone2(p=1, n=2)
          annotation (
            Placement(
              transformation(
                extent={{-160,-90},{-140,-70}})));

        BondGraph.BG_linear.Elements.Jone jone3(p=2, n=4)
          annotation (
            Placement(
              transformation(
                extent={{40,-60},{60,-40}})));

        BondGraph.BG_linear.Elements.Jone jone4(p=1, n=2)
          annotation (
            Placement(
              transformation(
                extent={{140,40},{160,60}})));

        BondGraph.BG_linear.Elements.Jzero jzero1(p=1, n=2)
          annotation (
            Placement(
              transformation(
                extent={{-120,70},{-100,90}})));

        BondGraph.BG_linear.Elements.Jzero jzero2(p=1, n=2)
          annotation (
            Placement(
              transformation(
                extent={{-120,-90},{-100,-70}})));

        BondGraph.BG_linear.Elements.Jzero jzero3(p=1, n=2)
          annotation (
            Placement(
              transformation(
                extent={{100,-60},{120,-40}})));

      equation
        connect(jone1.port_n[1], jzero1.port_p[1]) annotation (Line(
            points={{-141,79.5},{-129.5,79.5},{-129.5,80},{-119,80}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jone1.port_n[2], r_h_1.port) annotation (Line(
            points={{-141,80.5},{-141,50}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero1.port_n[2], c_h_1.port) annotation (Line(
            points={{-101,80.5},{-101,50}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(port_p1, jone1.port_p[1]) annotation (Line(
            points={{-190,80},{-159,80}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero1.port_n[1], tf1.port_p) annotation (Line(
            points={{-101,79.5},{-90.5,79.5},{-90.5,80},{-79,80}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jone2.port_n[1],jzero2. port_p[1]) annotation (Line(
            points={{-141,-80.5},{-130,-80.5},{-130,-80},{-119,-80}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jone2.port_n[2], r_h_2.port) annotation (Line(
            points={{-141,-79.5},{-141,-50}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero2.port_n[2], c_h_2.port) annotation (Line(
            points={{-101,-79.5},{-101,-50}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(port_p2,jone2. port_p[1]) annotation (Line(
            points={{-190,-80},{-159,-80}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero2.port_n[1], tf2.port_p) annotation (Line(
            points={{-101,-80.5},{-90.5,-80.5},{-90.5,-80},{-79,-80}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(tf1.port_n, sensor_i_f1.port_p) annotation (Line(
            points={{-61,80},{-50,80},{-50,-50},{-39,-50}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_i_f1.port_n, jone3.port_p[1]) annotation (Line(
            points={{-21,-50},{10,-50},{10,-50.5},{41,-50.5}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_i_f1.port_out, deadZone.u)  annotation (Line(
            points={{-30,-59},{-30,-70},{-10,-70},{-10,20},{-2,20}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(deadZone.y, mse1.port_in)  annotation (Line(
            points={{21,20},{50,20},{50,-1}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(jone3.port_n[2], i_m_1.port) annotation (Line(
            points={{59,-50.25},{80,-50.25},{80,80},{59,80}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jone3.port_n[3], r_m_1.port) annotation (Line(
            points={{59,-49.75},{70,-49.75},{70,50},{59,50}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jone3.port_n[4], mse1.port) annotation (Line(
            points={{59,-49.25},{59,-10}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jone3.port_n[1], jzero3.port_p[1]) annotation (Line(
            points={{59,-50.75},{80.5,-50.75},{80.5,-50},{101,-50}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero3.port_n[1], port_n) annotation (Line(
            points={{119,-50.5},{170,-50.5},{170,0},{190,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jone4.port_n[2],r_m_2. port) annotation (Line(
            points={{159,50.5},{159,80}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero3.port_n[2], jone4.port_p[1]) annotation (Line(
            points={{119,-49.5},{120,-49.5},{120,50},{141,50}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jone4.port_n[1], c_m_1.port) annotation (Line(
            points={{159,49.5},{159,20}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(tf2.port_n, jone3.port_p[2]) annotation (Line(
            points={{-61,-80},{30,-80},{30,-49.5},{41,-49.5}},
            color={0,0,0},
            smooth=Smooth.None));

        annotation (
          Diagram(
            coordinateSystem(
              preserveAspectRatio=true,
              extent={{-200,-100},{200,100}}),
            graphics),
          Icon(
            coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-200},{100,200}}),
            graphics={Rectangle(
                extent={{-100,200},{100,-140}},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Rectangle(
                extent={{-80,180},{80,-120}},
                lineColor={0,0,0},
                fillColor={125,190,255},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-80,32},{80,12}},
                lineColor={0,0,0},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-10,12},{10,-200}},
                lineColor={0,0,0},
                fillColor={175,175,175},
                fillPattern=FillPattern.Solid)}));
      end Cylinder_lin;
    end Cylinders;

    package Pumps
      model Pump_lin

        parameter Real eta_hm = 0.87 "hydraulic mechanical efficiency factor";
        parameter Real eta_vol = 1 "volumetric efficiency factor";
        parameter Real i_mech = 0.05
          "moment of inertia of moving pump parts reduced relative to rotational velocity of pump shaft";
        parameter Real r_mech = 1e-12
          "mechanical resistance of moving pump parts reduced relative to rotational velocity of pump shaft";
        parameter Real c_mech = 1e-6
          "mechanical compliance of pump parts reduced relative to rotational velocity of pump shaft";
        parameter Real v_disp = 16e-6
          "displaced volume of pump per shaft revolution";
        parameter Real r_leak = 1.25997480050399e12
          "inner pump leakage resistance";
        parameter Real f_mech_0 = 0 "initial rotational velosity of pump shaft";
        parameter Real e_mech_0 = 0 "initial torque of pump shaft";

        BondGraph.Interfaces.Port_p port_p
          annotation (
            Placement(
              transformation(
                extent={{-200,-10},{-180,10}}),
              iconTransformation(
                extent={{-100,-10},{-80,10}})));

        BondGraph.Interfaces.Port_n port_n
          annotation (
            Placement(
              transformation(
                extent={{180,-10},{200,10}}),
              iconTransformation(
                extent={{80,-10},{100,10}})));

        BondGraph.BG_linear.Sensors.Sensor_p_e sensor_p_e1(gain=1 - eta_hm)
          annotation (
            Placement(
              transformation(
                extent={{-160,-10},{-140,10}})));

        BondGraph.BG_linear.Sensors.Sensor_p_f sensor_p_f1(gain=1 - eta_vol)
          annotation (
            Placement(
              transformation(
                extent={{100,-10},{120,10}})));

        BondGraph.BG_linear.Elements.MSe mse1
          annotation (
            Placement(
              transformation(
                extent={{-120,-50},{-100,-30}})));

        BondGraph.BG_linear.Elements.MSf msf1
          annotation (
            Placement(
              transformation(
                extent={{140,-50},{160,-30}})));

        BondGraph.BG_linear.Elements.C c1(c=c_mech, e_0=e_mech_0)
          annotation (
            Placement(
              transformation(
                extent={{-40,30},{-20,50}})));

        BondGraph.BG_linear.Elements.I i1(i=i_mech, f_0=f_mech_0)
          annotation (
            Placement(
              transformation(
                extent={{-80,30},{-60,50}})));

        BondGraph.BG_linear.Elements.R r1(r=r_mech)
          annotation (
            Placement(
              transformation(
                extent={{-80,-50},{-60,-30}})));

        BondGraph.BG_linear.Elements.R r2(r=r_leak)
          annotation (
            Placement(
              transformation(
                extent={{60,-50},{80,-30}})));

        BondGraph.BG_linear.Elements.TF tf1(r=v_disp/(4*asin(1)))
          annotation (
            Placement(
              transformation(
                extent={{20,-10},{40,10}})));

        BondGraph.BG_linear.Elements.Jone jone1(p=1, n=2)
          annotation (
            Placement(
              transformation(
                extent={{-120,-10},{-100,10}})));

        BondGraph.BG_linear.Elements.Jone jone2(p=1, n=3)
          annotation (
            Placement(
              transformation(
                extent={{-80,-10},{-60,10}})));

        BondGraph.BG_linear.Elements.Jzero jzero1(n=2, p=1)
          annotation (
            Placement(
              transformation(
                extent={{-40,-10},{-20,10}})));

        BondGraph.BG_linear.Elements.Jzero jzero2(n=2, p=1)
          annotation (
            Placement(
              transformation(
                extent={{60,-10},{80,10}})));

        BondGraph.BG_linear.Elements.Jzero jzero3(n=2, p=1)
          annotation (
            Placement(
              transformation(
                extent={{140,-10},{160,10}})));

      equation
        connect(port_p, sensor_p_e1.port_p) annotation (Line(
            points={{-190,0},{-159,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_p_e1.port_n, jone1.port_p[1]) annotation (Line(
            points={{-141,0},{-119,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jone1.port_n[1], jone2.port_p[1]) annotation (Line(
            points={{-101,-0.5},{-89.5,-0.5},{-89.5,0},{-79,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jone2.port_n[1], jzero1.port_p[1]) annotation (Line(
            points={{-61,-0.666667},{-50.5,-0.666667},{-50.5,0},{-39,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero1.port_n[1], tf1.port_p) annotation (Line(
            points={{-21,-0.5},{0,-0.5},{0,0},{21,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(tf1.port_n, jzero2.port_p[1]) annotation (Line(
            points={{39,0},{61,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero2.port_n[1], sensor_p_f1.port_p) annotation (Line(
            points={{79,-0.5},{90.5,-0.5},{90.5,0},{101,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_p_f1.port_n, jzero3.port_p[1]) annotation (Line(
            points={{119,0},{141,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero3.port_n[1], port_n) annotation (Line(
            points={{159,-0.5},{190,-0.5},{190,0}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jone1.port_n[2], mse1.port) annotation (Line(
            points={{-101,0.5},{-101,-40}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jone2.port_n[2], i1.port) annotation (Line(
            points={{-61,5.55112e-017},{-61,40}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jone2.port_n[3], r1.port) annotation (Line(
            points={{-61,0.666667},{-61,-40}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero1.port_n[2], c1.port) annotation (Line(
            points={{-21,0.5},{-21,40}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero2.port_n[2], r2.port) annotation (Line(
            points={{79,0.5},{79,-40}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(jzero3.port_n[2], msf1.port) annotation (Line(
            points={{159,0.5},{159,-40}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_p_e1.port_out, mse1.port_in) annotation (Line(
            points={{-150,-9},{-150,-20},{-110,-20},{-110,-31}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(sensor_p_f1.port_out, msf1.port_in) annotation (Line(
            points={{110,-9},{110,-20},{150,-20},{150,-31}},
            color={0,0,0},
            smooth=Smooth.None));

        annotation (
          Diagram(
            coordinateSystem(
              preserveAspectRatio=true,
              extent={{-200,-100},{200,100}}),
            graphics),
          Icon(
            coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}}),
            graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={125,190,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}), Text(
                extent={{0,30},{0,-30}},
                pattern=LinePattern.None,
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                textString="pump
lin")}));
      end Pump_lin;
    end Pumps;
  end Examples;

  package Blocks

    block Const "constant matrix"

      BondGraph.Interfaces.Port_out port_out[size(c, 1),size(c, 2)]
        annotation (
          Placement(
            transformation(
              extent={{-10,-100},{10,-80}})));

      parameter Real c[:,:] = {{0}} "constant matrix";

    equation
      port_out = c;

      annotation (
        Icon(
          graphics={Ellipse(
              extent={{-100,100},{100,-100}},
              fillColor={51,249,65},
              fillPattern=FillPattern.Solid,
              lineColor={0,0,0}), Text(
              extent={{0,30},{0,-30}},
              pattern=LinePattern.None,
              fillColor={85,170,255},
              fillPattern=FillPattern.Solid,
              lineColor={0,0,0},
              textString="const")}));
    end Const;

    block Dis_der "dynamic element with discrete der(state)"

      import BondGraph.Functions.*;

      BondGraph.Interfaces.Port_in port_in[size(x_der, 1)]
        annotation (
          Placement(
            transformation(
              extent={{-10,80},{10,100}})));

      BondGraph.Interfaces.Port_out port_out[size(x_der, 1)]
        annotation (
          Placement(
            transformation(
              extent={{-10,-100},{10,-80}})));

      parameter Real x_der[:,2] = {{-1,1}} "discrete der(state) values";
      parameter Real x_0[size(x_der,1)] = {0} "initial state";
      parameter Real u_diff = 1e-3 "switching threshold";
      parameter Integer f_type = 1 "dynamic function type";

    protected
      Real u[size(x_der,1)] "input";
      Real x[size(x_der,1)] "state";
      Real t_u "switch time";
      Real u_u[size(x_der,1)] "switch input";
      Real x_u[size(x_der,1)] "switch state";
      Real x_der_u[size(x_der,1)] "switch der(state)";
      Real x_u_min[size(x_der,1)] "switch min state";
      Real x_u_max[size(x_der,1)] "switch max state";

    algorithm
      when {initial(),(u-u_u)*(u-u_u)>u_diff^2} then
        t_u := time;
        u_u := u;
        x_u := pre(x);
        if f_type==1 then
          x_der_u := 0.5
                     *
                     transpose(
                       x_der
                       *
                       { -sign(u_u-x_u)+ones(size(x_der,1)),
                         sign(u_u-x_u)+ones(size(x_der,1))})
                     *
                     {1};
        else
          x_der_u := 0.5
                     *
                     transpose(
                       x_der
                       *
                       { -sign(abs(u_u)-abs(x_u))+ones(size(x_der,1)),
                         sign(abs(u_u)-abs(x_u))+ones(size(x_der,1))})
                     *
                     {1};
        end if;
        for i in 1:size(x_der,1) loop
          x_u_min[i] := min({u_u[i],x_u[i]});
          x_u_max[i] := max({u_u[i],x_u[i]});
        end for;
      end when;

    initial equation
      pre(x) = x_0;

    equation
      port_in = u;
      port_out = x;

      x = saturation(x_der_u*(time-t_u)+x_u, x_u_min, x_u_max);

      annotation (
        Icon(
          graphics={Ellipse(
              extent={{-100,100},{100,-100}},
              fillColor={51,249,65},
              fillPattern=FillPattern.Solid,
              lineColor={0,0,0}), Text(
              extent={{0,30},{0,-30}},
              pattern=LinePattern.None,
              fillColor={85,170,255},
              fillPattern=FillPattern.Solid,
              lineColor={0,0,0},
              textString="dis
der")}));
    end Dis_der;

    block Dis_der_int "dynamic element with discrete der(state)"

      import BondGraph.Functions.*;

      BondGraph.Interfaces.Port_in port_in[size(x_der, 1)]
        annotation (
          Placement(
            transformation(
              extent={{-10,80},{10,100}})));

      BondGraph.Interfaces.Port_out port_out[size(x_der, 1)]
        annotation (
          Placement(
            transformation(
              extent={{-10,-100},{10,-80}})));

      parameter Real x_der[:,2] = {{-1,1}} "discrete der(state) values";
      parameter Real x_0[size(x_der,1)] = {0} "initial state";
      parameter Real u_diff = 1e-3 "switching threshold";
      parameter Integer f_type = 1 "dynamic function type";

    protected
      Real u[size(x_der,1)] "input";
      Real x[size(x_der,1)] "state";
      Real x_sat[size(x_der,1)] "limited state";
      Real u_u[size(x_der,1)] "switch input";
      Real x_u[size(x_der,1)] "switch state";
      Real x_der_u[size(x_der,1)] "switch der(state)";
      Real x_u_min[size(x_der,1)] "switch min state";
      Real x_u_max[size(x_der,1)] "switch max state";

    algorithm
      when {initial(),(u-u_u)*(u-u_u)>u_diff^2} then
        u_u := u;
        x_u := pre(x_sat);
        if f_type==1 then
          x_der_u := 0.5
                     *
                     transpose(
                       x_der
                       *
                       { -sign(u_u-x_u)+ones(size(x_der,1)),
                         sign(u_u-x_u)+ones(size(x_der,1))})
                     *
                     {1};
        else
          x_der_u := 0.5
                     *
                     transpose(
                       x_der
                       *
                       { -sign(abs(u_u)-abs(x_u))+ones(size(x_der,1)),
                         sign(abs(u_u)-abs(x_u))+ones(size(x_der,1))})
                     *
                     {1};
        end if;
        for i in 1:size(x_der,1) loop
          x_u_min[i] := min({u_u[i],x_u[i]});
          x_u_max[i] := max({u_u[i],x_u[i]});
        end for;
        reinit(x,x_u);
      end when;

    initial equation
      x = x_0;
      pre(x_sat) = x_0;

    equation
      port_in = u;
      port_out = x_sat;

      x_sat = saturation(x, x_u_min, x_u_max);
      der(x) = x_der_u;

      annotation (
        Icon(
          graphics={Ellipse(
              extent={{-100,100},{100,-100}},
              fillColor={51,249,65},
              fillPattern=FillPattern.Solid,
              lineColor={0,0,0}), Text(
              extent={{0,30},{0,-30}},
              pattern=LinePattern.None,
              fillColor={85,170,255},
              fillPattern=FillPattern.Solid,
              lineColor={0,0,0},
              textString="dis
der")}));
    end Dis_der_int;

    block Lin_sat "linear dynamic element with limited state and der(state)"

      import BondGraph.Functions.*;

      BondGraph.Interfaces.Port_in port_in[size(b, 2)]
        annotation (
          Placement(
            transformation(
              extent={{-10,80},{10,100}})));

      BondGraph.Interfaces.Port_out port_out[size(b, 1)]
        annotation (
          Placement(
            transformation(
              extent={{-10,-100},{10,-80}})));

      parameter Real b[:,:] = {{1}} "input gain matrix";
      parameter Real a[size(b,1),size(b,1)] = {{-1}} "state gain matrix";
      parameter Real x_0[size(b,1)] = {0} "initial state";
      parameter Real x_min[size(b,1)] = {-1} "minimal state values";
      parameter Real x_max[size(b,1)] = {1} "maximal state values";
      parameter Real x_der_min[size(b,1)] = {-1} "minimal der(state) values";
      parameter Real x_der_max[size(b,1)] = {1} "maximal der(state) values";

    protected
      Real u[size(b,2)] "input";
      Real x[size(b,1)](start = x_0) "state";

    equation
      port_in = u;
      port_out = saturation(x, x_min, x_max);

      der(x) = saturation(a*x + b*u, x_der_min, x_der_max);

      annotation (
        Icon(
          graphics={Ellipse(
              extent={{-100,100},{100,-100}},
              fillColor={51,249,65},
              fillPattern=FillPattern.Solid,
              lineColor={0,0,0}), Text(
              extent={{0,30},{0,-30}},
              pattern=LinePattern.None,
              fillColor={85,170,255},
              fillPattern=FillPattern.Solid,
              lineColor={0,0,0},
              textString="lin
sat")}));
    end Lin_sat;

    block Prop_sat "proportional element with limited state"

      import BondGraph.Functions.*;

      BondGraph.Interfaces.Port_in port_in[size(b, 2)]
        annotation (
          Placement(
            transformation(
              extent={{-10,80},{10,100}})));

      BondGraph.Interfaces.Port_out port_out[size(b, 1)]
        annotation (
          Placement(
            transformation(
              extent={{-10,-100},{10,-80}})));

      parameter Real b[:,:] = {{1}} "input gain matrix";
      parameter Real x_min[size(b,1)] = {-1} "minimal state values";
      parameter Real x_max[size(b,1)] = {1} "maximal state values";

    protected
      Real u[size(b,2)] "input";
      Real x[size(b,1)] "state";

    equation
      port_in = u;
      port_out = x;

      x = saturation(b*u, x_min, x_max);

      annotation (
        Icon(
          graphics={Ellipse(
              extent={{-100,100},{100,-100}},
              fillColor={51,249,65},
              fillPattern=FillPattern.Solid,
              lineColor={0,0,0}), Text(
              extent={{0,30},{0,-30}},
              pattern=LinePattern.None,
              fillColor={85,170,255},
              fillPattern=FillPattern.Solid,
              lineColor={0,0,0},
              textString="prop
sat")}));
    end Prop_sat;

    block Sig_c "control signal"

      import BondGraph.Functions.*;

      BondGraph.Interfaces.Port_in port_in
        annotation (
          Placement(
            transformation(
              extent={{-10,80},{10,100}})));

      BondGraph.Interfaces.Port_out port_out[size(p_x, 1)]
        annotation (
          Placement(
            transformation(
              extent={{-10,-100},{10,-80}})));

      parameter Real p_x[:] = {0,1} "states assignment";
      parameter Real p_cn[size(p_x,1)] = {1,1}
        "coverage parameter for state < p_x[i]";
      parameter Real p_cp[size(p_x,1)] = {1,1}
        "coverage parameter for state > p_x[i]";

    protected
      parameter Integer m = size(p_x,1);
      parameter Real p_xn[m] = cat(1, {1}, ones(m-1) ./ (p_x[2:m] - p_x[1:m-1]));
      parameter Real p_xp[m] = cat(1, ones(m-1) ./ (p_x[2:m] - p_x[1:m-1]), {1});

      Real x "position";
      Real s_c[m] "control signal";

    equation
      port_in = x;
      port_out = s_c;

      s_c = ones(m)
            - saturation(-p_xn ./ p_cn .* (x * ones(m) - p_x), zeros(m), ones(m))
            - saturation(+p_xp ./ p_cp .* (x * ones(m) - p_x), zeros(m), ones(m));

      annotation (
        Icon(
          graphics={Ellipse(
              extent={{-100,100},{100,-100}},
              fillColor={51,249,65},
              fillPattern=FillPattern.Solid,
              lineColor={0,0,0}), Text(
              extent={{0,30},{0,-30}},
              pattern=LinePattern.None,
              fillColor={85,170,255},
              fillPattern=FillPattern.Solid,
              lineColor={0,0,0},
              textString="sig
c")}));
    end Sig_c;

    block Sig_cs "control signal of specialized characteristic"

      BondGraph.Interfaces.Port_in port_in
        annotation (
          Placement(
            transformation(
              extent={{-10,80},{10,100}})));

      BondGraph.Interfaces.Port_out port_out
        annotation (
          Placement(
            transformation(
              extent={{-10,-100},{10,-80}})));

      parameter Real w[4] = {1,0,0,0}
        "weights for {polynomial, exponential, hyperbolic, s-type} characteristic, w in R^4 and w[i] in [0, +inf)";
      parameter Real p[4] = {1,1,1,1}
        "parameters for {polynomial, exponential, hyperbolic, s-type} characteristic, p in R^4 and p[i] in (0, +inf)";
      parameter Real s_cs_min = 1e-9
        "minimal control signal value, s_cs_min in R and s_cs_min in (0, 1)";
      parameter Real s_cs_max = 1
        "maximal control signal value, s_cs_max in R and s_cs_max in (s_cs_min, 1]";

    protected
      Real s_c "control signal, s_c in R and s_c in [0,1]";
      Real s_cs
        "control signal of specialized characteristic, s_cs in R and s_cs in [s_cs_min, s_cs_max]";

    equation
      port_in = s_c;
      port_out = s_cs;

      s_cs = (w / (ones(4) * w))
             *
             { (s_cs_max - s_cs_min) * abs(s_c) ^ p[1] + s_cs_min,
               exp(log(s_cs_max / s_cs_min) * abs(s_c) ^ p[2] + log(s_cs_min)),
               (s_cs_max - s_cs_min) * (2/(2 - abs(s_c) ^ p[3]) - 1) + s_cs_min,
               (s_cs_max-s_cs_min) * 0.5 * (tanh(2 * p[4] * s_c - p[4]) / tanh(p[4]) + 1) + s_cs_min};

      annotation (
        Icon(
          graphics={Ellipse(
              extent={{-100,100},{100,-100}},
              fillColor={51,249,65},
              fillPattern=FillPattern.Solid,
              lineColor={0,0,0}), Text(
              extent={{0,30},{0,-30}},
              pattern=LinePattern.None,
              fillColor={85,170,255},
              fillPattern=FillPattern.Solid,
              lineColor={0,0,0},
              textString="sig
cs
")}));
    end Sig_cs;
  end Blocks;

  package Functions

    function saturation

      input Real x;
      input Real x_min;
      input Real x_max;
      output Real x_sat;

    algorithm
      if x < x_min then
        x_sat := x_min;
      elseif x > x_max then
        x_sat := x_max;
      else
        x_sat := x;
      end if;

      annotation (
        derivative(order=1)=saturation_d1);
    end saturation;

    function saturation_d1

      input Real x;
      input Real x_min;
      input Real x_max;
      input Real x_d1;
      input Real x_min_d1;
      input Real x_max_d1;
      output Real x_sat_d1;

    algorithm
      if x < x_min then
        x_sat_d1 := 1 * x_min_d1;
      elseif x > x_max then
        x_sat_d1 := 1 * x_max_d1;
      else
        x_sat_d1 := 1 * x_d1;
      end if;

      annotation (
        derivative(order=1)=saturation_d2);
    end saturation_d1;

    function saturation_d2

      input Real x;
      input Real x_min;
      input Real x_max;
      input Real x_d1;
      input Real x_min_d1;
      input Real x_max_d1;
      input Real x_d2;
      input Real x_min_d2;
      input Real x_max_d2;
      output Real x_sat_d2;

    algorithm
      x_sat_d2 := 0;

      annotation (
        derivative(order=1)=saturation_d3);
    end saturation_d2;

    function saturation_d3

      input Real x;
      input Real x_min;
      input Real x_max;
      input Real x_d1;
      input Real x_min_d1;
      input Real x_max_d1;
      input Real x_d2;
      input Real x_min_d2;
      input Real x_max_d2;
      input Real x_d3;
      input Real x_min_d3;
      input Real x_max_d3;
      output Real x_sat_d3;

    algorithm
      x_sat_d3 := 0;

    end saturation_d3;
  end Functions;

  package Interfaces

    connector Port
      Real e;
      Real f;

      annotation (
        Icon(
          graphics={Ellipse(
              extent={{-100,100},{100,-100}},
              fillColor={135,135,135},
              fillPattern=FillPattern.Solid,
              lineColor={0,0,0})}));
    end Port;

    connector Port_in = input Real
      annotation (
        Icon(
          graphics={Polygon(
            points={{-100,100},{0,-100},{100,100},{-100,100}},
            lineColor={0,0,0},
            smooth=Smooth.None,
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid)}));
    connector Port_n
      Real e;
      Real f;

      annotation (
        Icon(
          graphics={Ellipse(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={255,0,0},
              fillPattern=FillPattern.Solid)}));
    end Port_n;

    connector Port_out = output Real "output"
      annotation (
        Icon(
          graphics={Polygon(
            points={{-100,100},{0,-100},{100,100},{-100,100}},
            lineColor={0,0,0},
            smooth=Smooth.None,
            fillColor={255,0,0},
            fillPattern=FillPattern.Solid)}));
    connector Port_p
      Real e;
      Real f;

      annotation (
        Icon(
          graphics={Ellipse(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={0,0,255},
              fillPattern=FillPattern.Solid)}));
    end Port_p;
  end Interfaces;
  annotation (
    Icon(
      graphics),                            uses(Modelica(version="3.2")));
end BondGraph;
